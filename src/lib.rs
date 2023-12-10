use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::mem::swap;
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
    x: usize,
}

impl LCG {
    const A: usize = 1103515245;
    const C: usize = 12345;

    pub fn new(seed: usize) -> Self {
        LCG { x: seed }
    }

    pub fn gen(&mut self) -> usize {
        self.x = Self::A.wrapping_mul(self.x).wrapping_add(Self::C);
        self.x
    }

    pub fn gen_by(&mut self, r: Range<usize>) -> usize {
        let n = r.end - r.start + 1;
        r.start + self.gen() % n // ignore modulo bias
    }
}

pub fn inc_mod(i: usize, n: usize) -> usize {
    let j = i + 1;
    if j >= n {
        j - n
    } else {
        j
    }
}

pub fn dec_mod(i: usize, n: usize) -> usize {
    if i > 0 {
        i - 1
    } else {
        n - 1
    }
}

#[derive(Clone)]
pub struct TourData {
    seq: Vec<NodeId>,
    pos: Vec<usize>,
}

impl TourData {
    pub fn new(seq: Vec<NodeId>) -> Self {
        let n = seq.len();
        let mut x = Self {
            seq,
            pos: vec![0; n],
        };
        x.repos();
        x
    }

    pub fn repos(&mut self) {
        let Self { seq, pos } = self;
        for (i, &x) in seq.iter().enumerate() {
            pos[x] = i;
        }
    }

    pub fn n(&self) -> usize {
        self.seq.len()
    }

    pub fn inc(&self, i: usize) -> usize {
        let n = self.n();
        inc_mod(i, n)
    }

    pub fn dec(&self, i: usize) -> usize {
        let n = self.n();
        dec_mod(i, n)
    }

    pub fn next_id(&self, t: NodeId) -> NodeId {
        self.seq[self.inc(self.pos[t])]
    }

    pub fn prev_id(&self, t: NodeId) -> NodeId {
        self.seq[self.dec(self.pos[t])]
    }

    pub fn id(&self, t: NodeId, dir: bool) -> NodeId {
        if dir {
            self.next_id(t)
        } else {
            self.prev_id(t)
        }
    }

    pub fn seglen(&self, t1: NodeId, t2: NodeId) -> usize {
        let i1 = self.pos[t1];
        let i2 = self.pos[t2];
        let n = self.n();
        if i1 <= i2 {
            i2 - i1
        } else {
            i2 + n - i1
        }
    }

    pub fn reverse(&mut self, x: NodeId, y: NodeId) {
        let m = (self.seglen(x, y) + 1) / 2;
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

    pub fn opt2move(&mut self, t1: NodeId, t2: NodeId, t3: NodeId, t4: NodeId) {
        let n = self.n();
        if self.seglen(t2, t3) < n / 2 {
            self.reverse(t2, t3);
        } else {
            self.reverse(t4, t1);
        }
    }

    pub fn perturb(&mut self, rng: &mut LCG, buf: &mut Self) -> Option<[NodeId; 8]> {
        let n = self.n();
        let Self { seq, .. } = self;

        if n > 8 {
            let i1 = rng.gen_by(1..n / 4);
            let i2 = inc_mod(i1, n);
            let i3 = rng.gen_by(1..n / 4) + i1;
            let i4 = inc_mod(i3, n);
            let i5 = rng.gen_by(1..n / 4) + i3;
            let i6 = inc_mod(i5, n);
            let i7 = n - 1;
            let i8 = 0;
            let ans = [
                seq[i1], seq[i2], seq[i3], seq[i4], seq[i5], seq[i6], seq[i7], seq[i8],
            ];

            buf.seq.clear();
            for slice in [&seq[..i1], &seq[i5..], &seq[i3..i5], &seq[i1..i3]] {
                buf.seq.extend_from_slice(slice);
            }
            buf.repos();
            swap(self, buf);

            Some(ans)
        } else {
            None
        }
    }
}

#[derive(Clone)]
struct Optim {
    q: VecDeque<NodeId>,
    dlb: Vec<bool>,
}

impl Optim {
    pub fn new(n: usize) -> Self {
        Self {
            q: (0..n).collect(),
            dlb: vec![false; n],
        }
    }

    fn reset_only(&mut self, x: NodeId) {
        let Self { q, dlb } = self;
        if dlb[x] {
            dlb[x] = false;
            q.push_back(x);
        }
    }

    pub fn reset(&mut self, data: &TourData, s: &[NodeId]) {
        for &x in s {
            self.reset_only(x);
            self.reset_only(data.next_id(x));
            self.reset_only(data.prev_id(x));
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum StepError {
    EmptyQueue,
    HasDLB,
    NoGainMove,
}

#[derive(Clone)]
pub struct Tour<'a> {
    cost: i64,
    tsp: &'a TSP,
    data: TourData,
    opt2: Optim,
    opt3: Optim,
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
                        .min_by_key(|&i| tsp.d_id(a, i))
                })
                .unwrap_or(0);
            seq.push(x);
            used[x] = true;
        }

        Self {
            cost: tsp.compute_cost(&seq),
            tsp,
            data: TourData::new(seq),
            opt2: Optim::new(tsp.n),
            opt3: Optim::new(tsp.n),
        }
    }

    pub fn seq(&self) -> &[NodeId] {
        &self.data.seq
    }

    pub fn cost(&self) -> i64 {
        self.cost
    }

    pub fn perturb(&mut self, rng: &mut LCG, buf: &mut TourData) {
        let Self {
            tsp, data, cost, opt2, ..
        } = self;
        if let Some(indices) = data.perturb(rng, buf) {
            opt2.reset(data, &indices);
            *cost = tsp.compute_cost(&data.seq);
        }
    }

    pub fn opt2step(&mut self) -> Result<i64, StepError> {
        let Self {
            tsp,
            data,
            opt2,
            cost,
            ..
        } = self;
        let Optim { q, dlb } = opt2;
        let d = |i, j| tsp.d_id(i, j);
        let t1 = q.pop_front().ok_or(StepError::EmptyQueue)?;
        if dlb[t1] {
            Err(StepError::HasDLB)
        } else {
            dlb[t1] = true;
            for dir in [true, false] {
                let t2 = data.id(t1, dir);
                for &t3 in tsp.c[t1].iter().filter(|&&x| !dlb[x]) {
                    if d(t1, t2) < d(t1, t3) {
                        break;
                    }
                    let t4 = data.id(t3, dir);

                    if t1 != t4 && t2 != t3 {
                        let g = d(t1, t2) - d(t1, t3) + d(t3, t4) - d(t2, t4);
                        if g > 0 {
                            if dir {
                                data.opt2move(t1, t2, t3, t4);
                            } else {
                                data.opt2move(t2, t1, t4, t3);
                            }
                            opt2.reset(data, &[t1, t2, t3, t4]);
                            *cost -= g;
                            return Ok(g);
                        }
                    }
                }
            }
            Err(StepError::NoGainMove)
        }
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
            duration: duration - Duration::from_millis(200),
        }
    }

    pub fn ok(&self) -> bool {
        self.t.elapsed() < self.duration
    }
}

pub struct Solver<'a> {
    timer: Timer,
    rng: LCG,
    cur: Tour<'a>,
    bst: Tour<'a>,
}

impl<'a> Solver<'a> {
    pub fn new(tsp: &'a TSP, duration: Duration, seed: usize) -> Self {
        let t = Tour::new(&tsp);
        Self {
            timer: Timer::new(duration),
            rng: LCG::new(seed),
            cur: t.clone(),
            bst: t,
        }
    }

    pub fn run(&mut self) -> usize {
        let Self {
            timer,
            cur,
            rng,
            bst,
        } = self;
        let n = cur.tsp.n;
        let mut buf = TourData {
            seq: Vec::with_capacity(n),
            pos: vec![0; n],
        };

        let mut cnt = 0;
        while timer.ok() {
            cnt += 1;
            cur.perturb(rng, &mut buf);

            while timer.ok() {
                if cur.opt2step() == Err(StepError::EmptyQueue) {
                    break;
                }
            }

            // Opt3::new(&mut self.cur)
            //     .take_while(|_| self.timer.ok());
            if cur.cost < bst.cost {
                *bst = cur.clone();
            } else {
                *cur = bst.clone();
            }
        }
        cnt
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
            .flat_map(|s| p.iter().map(move |t| d_pt(s, t).round() as i64))
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

    pub fn d_id(&self, i: usize, j: usize) -> i64 {
        self.m[i * self.n + j]
    }

    pub fn compute_cost(&self, seq: &[NodeId]) -> i64 {
        let mut cost = 0;
        for (i, &x) in seq.iter().enumerate() {
            let y = *seq.get(i + 1).unwrap_or(&seq[0]);
            cost += self.d_id(x, y);
        }
        cost
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

        let a = (999999.0, -999999.0);
        let b = (-999999.0, 999999.0);
        assert_eq!(d_pt(&a, &b) as i64, 2828424);
    }

    #[test]
    fn test_tsp_new() {
        let tsp = TSP::new(&POINTS);
        for i in 0..10 {
            let j = (i + 1) % 10;
            let x = &POINTS[i];
            let y = &POINTS[j];
            assert_eq!(tsp.d_id(i, j), d_pt(x, y).round() as i64);
            assert_eq!(tsp.d_id(j, i), d_pt(y, x).round() as i64);
        }
    }

    #[test]
    fn test_tsp_parse_tsp() {
        let path = Path::new("src/data/sample.tsp");
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
        assert_eq!(t.data.seq, SEQ);
    }

    #[test]
    fn test_opt2() {
        let tsp = TSP::new(&POINTS);
        let mut t = Tour::new(&tsp);
        for _ in 0..10 {
            let _ = t.opt2step();
        }
        assert_eq!(t.cost(), 276);
        assert_eq!(t.seq(), [0, 8, 4, 5, 3, 1, 7, 9, 6, 2]);
    }
}
