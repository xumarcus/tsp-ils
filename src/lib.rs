use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::ops::Range;
use std::time::{Duration, Instant};
type Point = (f64, f64);
type NodeId = usize;

pub fn d_pt(s: &Point, t: &Point) -> f64 {
    let dx = s.0 - t.0;
    let dy = s.1 - t.1;
    (dx * dx + dy * dy).sqrt()
}

pub struct LCG {
    x: u64,
}

impl LCG {
    const A: u64 = 1103515245;
    const C: u64 = 12345;

    pub fn new(seed: u64) -> Self {
        LCG { x: seed }
    }

    pub fn gen(&mut self) -> u64 {
        self.x = Self::A * self.x + Self::C;
        self.x
    }
}

#[derive(Clone)]
pub struct Tour<'a> {
    tsp: &'a TSP,
    seq: Vec<NodeId>,
    pos: Vec<usize>,
    cost: i64,
}

impl<'a> Tour<'a> {
    pub fn new(tsp: &'a TSP) -> Self {
        let mut seq: Vec<NodeId> = vec![];
        let mut used = vec![false; tsp.n];
        while seq.len() < tsp.n {
            let x = seq
                .last()
                .and_then(|&a| {
                    (0..tsp.n)
                        .filter(|&i| !used[i] && i != a)
                        .min_by_key(|&i| tsp.d_ix(a, i))
                })
                .unwrap_or(0);
            seq.push(x);
            used[x] = true;
        }
        Self::from_seq(tsp, seq)
    }

    pub fn from_seq(tsp: &'a TSP, seq: Vec<NodeId>) -> Self {
        let mut pos: Vec<usize> = vec![0; tsp.n];
        let mut cost = 0;
        for (i, &x) in seq.iter().enumerate() {
            let y = *seq.get(i + 1).unwrap_or(&seq[0]);
            pos[x] = i;
            cost += tsp.d_ix(x, y);
        }
        Self {
            tsp,
            seq,
            pos,
            cost,
        }
    }

    pub fn inc(&self, i: usize) -> usize {
        let n = self.tsp.n;
        let j = i + 1;
        if j >= n {
            j - n
        } else {
            j
        }
    }

    pub fn dec(&self, i: usize) -> usize {
        let n = self.tsp.n;
        if i > 0 {
            i - 1
        } else {
            n - 1
        }
    }

    pub fn next_of(&self, t: NodeId) -> NodeId {
        self.seq[self.inc(self.pos[t])]
    }

    pub fn prev_of(&self, t: NodeId) -> NodeId {
        self.seq[self.dec(self.pos[t])]
    }

    pub fn seglen(&self, t1: NodeId, t2: NodeId) -> usize {
        let i1 = self.pos[t1];
        let i2 = self.pos[t2];
        let n = self.tsp.n;
        if i1 <= i2 {
            i2 - i1
        } else {
            i2 + n - i1
        }
    }

    pub fn reverse(&mut self, x: NodeId, y: NodeId) {
        let m = self.seglen(x, y).div_ceil(2);
        let mut i = self.pos[x];
        let mut j = self.pos[y];
        for _ in 0..m {
            self.seq.swap(i, j);
            self.pos[self.seq[i]] = i;
            self.pos[self.seq[j]] = j;
            i = self.inc(i);
            j = self.dec(j);
        }
    }

    pub fn seq(&self) -> &[NodeId] {
        &self.seq
    }

    pub fn cost(&self) -> i64 {
        self.cost
    }
}

struct Timer {
    t: Instant,
    duration: Duration,
}

impl Timer {
    pub fn new(duration: Duration) -> Self {
        Self {
            t: Instant::now(),
            duration,
        }
    }

    pub fn ok(&self) -> bool {
        self.t.elapsed() < self.duration
    }
}

struct Opt2<'a, 'b> {
    tour: &'b mut Tour<'a>,
    r: Range<NodeId>,
}

impl<'a, 'b> Opt2<'a, 'b> {
    pub fn new(tour: &'b mut Tour<'a>) -> Self {
        let n = tour.tsp.n;
        Self { tour, r: 0..n }
    }
}

impl<'a, 'b> Iterator for Opt2<'a, 'b> {
    type Item = i64;
    fn next(&mut self) -> Option<Self::Item> {
        let tour = &mut self.tour;
        let d = |i, j| tour.tsp.d_ix(i, j);
        let n = tour.tsp.n;

        while let Some(t1) = self.r.next() {
            for dir in [true, false] {
                let t2 = if dir {
                    tour.next_of(t1)
                } else {
                    tour.prev_of(t1)
                };
                for &t3 in tour.tsp.c[t1].iter() {
                    let g1 = d(t1, t2) - d(t1, t3);
                    if g1 < 0 {
                        break;
                    }

                    let t4 = if dir {
                        tour.next_of(t3)
                    } else {
                        tour.prev_of(t3)
                    };
                    if t1 != t4 && t2 != t3 {
                        let g2 = d(t3, t4) - d(t2, t4);
                        let g = g1 + g2;
                        if g > 0 {
                            if dir {
                                if tour.seglen(t2, t3) < n / 2 {
                                    tour.reverse(t2, t3);
                                } else {
                                    tour.reverse(t4, t1);
                                }
                            } else {
                                if tour.seglen(t1, t4) < n / 2 {
                                    tour.reverse(t1, t4);
                                } else {
                                    tour.reverse(t3, t2);
                                }
                            }
                            tour.cost -= g;
                            return Some(g);
                        }
                    }
                }
            }
        }
        None
    }
}

pub struct Solver<'a> {
    timer: Timer,
    // rng: LCG,
    cur: Tour<'a>,
    bst: Tour<'a>,
}

impl<'a> Solver<'a> {
    pub fn new(tsp: &'a TSP) -> Self {
        let t = Tour::new(&tsp);
        Self {
            timer: Timer::new(Duration::from_secs(2)),
            // rng: LCG::new(4234),
            cur: t.clone(),
            bst: t,
        }
    }

    pub fn run(&mut self) {
        while self.timer.ok() {
            // self.cur.perturb(&mut self.rng);
            Opt2::new(&mut self.cur)
                .take_while(|_| self.timer.ok())
                .for_each(drop);
            // Opt3::new(&mut self.cur)
            //     .take_while(|_| self.timer.ok());
            if self.cur.cost < self.bst.cost {
                self.bst = self.cur.clone();
            } else {
                self.cur = self.bst.clone();
            }
        }
    }

    pub fn solution(&self) -> &Tour {
        &self.bst
    }
}

pub struct TSP {
    n: usize,
    m: Vec<i64>,
    c: Vec<Vec<NodeId>>,
}

impl TSP {
    const NEIGHBOR_LIMIT: usize = 30;

    pub fn new(p: &[Point]) -> Self {
        let n = p.len();
        let m: Vec<i64> = p
            .iter()
            .flat_map(|s| p.iter().map(|t| d_pt(s, t).round() as i64))
            .collect();
        let c = (0..n)
            .map(|i| {
                let mut v: Vec<NodeId> = (0..n).filter(|&j| i != j).collect();
                v.sort_by_key(|&j| m[i * n + j]);
                v.truncate(Self::NEIGHBOR_LIMIT);
                v
            })
            .collect();
        Self { n, m, c }
    }

    fn parse_tsp(r: impl BufRead, is_kattis: bool) -> io::Result<Self> {
        let v: Vec<Option<Point>> = r
            .lines()
            .map(|line| {
                line.map(|s| {
                    let mut r = s
                        .split(" ")
                        .skip(if is_kattis { 0 } else { 1 })
                        .map(|x| x.parse::<f64>());
                    Some((r.next()?.ok()?, r.next()?.ok()?))
                })
            })
            .collect::<io::Result<_>>()?;
        let p = v.into_iter().flatten().collect::<Vec<_>>();
        Ok(Self::new(&p))
    }

    pub fn parse_tsp_file(f: File) -> io::Result<Self> {
        Self::parse_tsp(BufReader::new(f), false)
    }

    pub fn parse_kattis() -> io::Result<Self> {
        let stdin = io::stdin();
        Self::parse_tsp(stdin.lock(), true)
    }

    pub fn d_ix(&self, i: usize, j: usize) -> i64 {
        self.m[i * self.n + j]
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::*;

    const POINTS: [Point; 10] = [
        (95.0129, 61.5432),
        (23.1139, 79.1937),
        (60.6843, 92.1813),
        (48.5982, 73.8207),
        (89.1299, 17.6266),
        (76.2097, 40.5706),
        (45.6468, 93.5470),
        (1.8504, 91.6904),
        (82.1407, 41.0270),
        (44.4703, 89.3650),
    ];
    const SEQ: [usize; 10] = [0, 8, 5, 4, 3, 9, 6, 2, 1, 7];

    #[test]
    fn test_d_pt() {
        let ans = [24, 6, 26, 69, 16, 4, 15, 40, 25, 98];
        for (i, d) in ans.into_iter().enumerate() {
            let x = &POINTS[SEQ[i]];
            let y = &POINTS[SEQ[(i + 1) % 10]];
            assert_eq!(d_pt(x, y).round() as i64, d);
        }
    }

    #[test]
    fn test_tsp_new() {
        let tsp = TSP::new(&POINTS);
        for i in 0..10 {
            let j = (i + 1) % 10;
            let x = &POINTS[i];
            let y = &POINTS[j];
            assert_eq!(tsp.d_ix(i, j), d_pt(x, y).round() as i64);
            assert_eq!(tsp.d_ix(j, i), d_pt(y, x).round() as i64);
        }
    }

    #[test]
    fn test_tsp_parse_tsp() {
        let path = Path::new("tests/sample.tsp");
        let tsp = File::open(path).and_then(TSP::parse_tsp_file).unwrap();
        let ans = [
            0, 74, 46, 48, 44, 28, 59, 98, 24, 58, 74, 0, 40, 26, 90, 66, 27, 25, 70, 24, 46, 40,
            0, 22, 80, 54, 15, 59, 55, 16, 48, 26, 22, 0, 69, 43, 20, 50, 47, 16, 44, 90, 80, 69,
            0, 26, 87, 114, 24, 85, 28, 66, 54, 43, 26, 0, 61, 90, 6, 58, 59, 27, 15, 20, 87, 61,
            0, 44, 64, 4, 98, 25, 59, 50, 114, 90, 44, 0, 95, 43, 24, 70, 55, 47, 24, 6, 64, 95, 0,
            61, 58, 24, 16, 16, 85, 58, 4, 43, 61, 0,
        ];
        assert_eq!(&tsp.m[..], ans);
    }

    #[test]
    fn test_tour_new() {
        let tsp = TSP::new(&POINTS);
        let t = Tour::new(&tsp);
        assert_eq!(t.seq, SEQ);
    }

    #[test]
    fn test_opt2() {
        let tsp = TSP::new(&POINTS);
        let mut t = Tour::new(&tsp);
        Opt2::new(&mut t).for_each(drop);
        assert_eq!(t.cost(), 276);
        assert_eq!(t.seq(), [0, 8, 4, 5, 3, 1, 7, 9, 6, 2]);
    }
}
