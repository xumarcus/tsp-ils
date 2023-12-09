use std::{fs::File, path::Path, time::Duration};

use tsp_ils::{Solver, TSP};

#[test]
fn main() {
    let t = File::open(Path::new("tests/data/berlin52.tsp"))
        .and_then(TSP::parse_tsp_file)
        .unwrap();
    let mut s = Solver::new(&t, Duration::from_secs(2), 4234);
    s.run();
    assert_eq!(s.solution().cost(), 7542);
}
