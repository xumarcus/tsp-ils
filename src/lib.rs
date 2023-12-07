use std::time::{Duration, Instant};
use std::path::Path;
use std::ops::Range;
use std::fs::File;
use std::io::{self, BufReader, BufRead};
type Point = (f64, f64);
type NodeId = usize;

pub fn d_pt(s: &Point, t: &Point) -> f64 {
    let dx = s.0 - t.0;
    let dy = s.1 - t.1;
    (dx * dx + dy * dy).sqrt()
}

pub struct LCG {
    x: u64
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
    cost: u64,
}

impl<'a> Tour<'a> {
    pub fn new(tsp: &'a TSP) -> Self {
        let mut seq: Vec<NodeId> = vec![];
        let mut used = vec![false; tsp.n];
        while seq.len() < tsp.n {
            let x = if let Some(&a) = seq.last() {
                (0..tsp.n)
                .filter(|&i| !used[i] && i != a)
                .min_by_key(|&i| tsp.d_ix(a, i))
                .unwrap()
            } else {
                0
            };
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
            println!("{} {} {}", tsp.d_ix(x, y), x, y);
        }
        Self { tsp, seq, pos, cost }
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

    pub fn reverse(&mut self, t1: NodeId, t2: NodeId) {
        let m = self.seglen(t1, t2).div_ceil(2);
        let mut i1 = self.pos[t1];
        let mut i2 = self.pos[t2];
        for _ in 0..m {
            self.seq.swap(i1, i2);
            self.pos[self.seq[i1]] = i1;
            self.pos[self.seq[i2]] = i2;
            i1 = self.inc(i1);
            i2 = self.dec(i2);
        }
    }

    pub fn cost(&self) -> u64 {
        self.cost
    }
}

struct Timer {
    t: Instant,
    duration: Duration,
}

impl Timer {
    pub fn new(duration: Duration) -> Self {
        Self { t: Instant::now(), duration }
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
    type Item = u64;
    fn next(&mut self) -> Option<Self::Item> {
        let tour = &mut self.tour;
        let d = |i, j| tour.tsp.d_ix(i, j);
        let n = tour.tsp.n;

        while let Some(t1) = self.r.next() {
            let t2 = tour.next_of(t1);
            for &t3 in tour.tsp.c[t1].iter() {
                let g1 = d(t1, t2) - d(t1, t3); 
                if g1 > 0 {
                    let t4 = tour.next_of(t3);
                    if t1 != t4 && t2 != t3 {
                        let g2 = g1 + d(t3, t4) - d(t2, t4);
                        if g2 > 0 {
                            if tour.seglen(t2, t3) < n / 2 {
                                tour.reverse(t2, t3);
                            } else {
                                tour.reverse(t4, t1);
                            }
                            tour.cost -= g2;
                            return Some(g2);
                        }
                    }
                } else {
                    break;  // g1 < 0 for all t3' after t3, as d(t1, t3') > d(t1, t3)
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
            let _ = Opt2::new(&mut self.cur)
                .take_while(|_| self.timer.ok());
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
    m: Vec<u64>,
    c: Vec<Vec<NodeId>>,
}

impl TSP {
    const NEIGHBOR_LIMIT: usize = 30;

    pub fn new(p: &[Point]) -> Self {
        let n = p.len();
        let m: Vec<u64> = p.iter() 
            .flat_map(|s| p.iter()
                .map(|t| d_pt(s, t) as u64))
            .collect();
        let c = (0..n)
            .map(|i| {
                let mut v: Vec<NodeId> = (0..n).collect();
                v.sort_by_key(|j| m[i * n + j]);
                v.truncate(Self::NEIGHBOR_LIMIT);
                v
            })
            .collect();
        Self {n, m, c}
    }

    fn pff(f: File) -> io::Result<Vec<Point>> {
        let v: Vec<Option<Point>> = BufReader::new(f)
            .lines()
            .map(|line| line
                .map(|s| {
                    let mut r = s.split(" ")
                        .skip(1)
                        .map(|x| x.parse::<f64>());
                    Some((r.next()?.ok()?, r.next()?.ok()?))
                })
            )
            .collect::<io::Result<_>>()?;
        Ok(v.into_iter().flatten().collect())
    }

    pub fn from_path(path: &Path) -> io::Result<Self> {
        let f = File::open(path)?;
        let p = Self::pff(f)?;
        Ok(Self::new(&p))
    }

    pub fn d_ix(&self, i: usize, j: usize) -> u64 {
        self.m[i * self.n + j]
    }
}

#[cfg(test)]
mod tests {
    use crate::*;

    const POINTS: [Point; 10] = [
        (95.0129, 61.5432),
        (23.1139, 79.1937),
        (60.6843, 92.1813),
        (48.5982, 73.8207),
        (89.1299, 17.6266),
        (76.2097, 40.5706),
        (45.6468, 93.5470),
        (1.8504 , 91.6904),
        (82.1407, 41.0270),
        (44.4703, 89.3650),
    ];
    const SEQ: [usize; 10] = [0, 8, 5, 4, 3, 9, 6, 2, 1, 7];

    #[test]
    fn test_d_pt() {
        let ds = [24, 5, 26, 69, 16, 4, 15, 39, 24, 97];
        for (i, d) in ds.into_iter().enumerate() {
            let x = &POINTS[SEQ[i]];
            let y = &POINTS[SEQ[(i + 1) % 10]];
            assert_eq!(d_pt(x, y) as u64, d);
        }
    }

    #[test]
    fn test_pff() {
        let p = File::open(Path::new("tests/sample.tsp"))
            .and_then(TSP::pff).unwrap();
        assert_eq!(&POINTS[..], &p[..]);
    }

    #[test]
    fn test_tour_new() {
        let tsp = TSP::new(&POINTS);
        let t = Tour::new(&tsp);
        assert_eq!(t.seq, SEQ);
    }
}