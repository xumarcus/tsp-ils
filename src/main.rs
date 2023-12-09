#[allow(unused_imports)]
use std::{io, fs::File, path::Path, time::Duration};

use tsp_ils::{Solver, TSP};

fn main() -> io::Result<()> {
    // let t = File::open(Path::new("src/data/berlin52.tsp")).and_then(TSP::parse_tsp_file)?;
    let t = TSP::parse_kattis()?;
    let mut s = Solver::new(&t, Duration::from_secs(2), 4234);
    s.run();
    for x in s.solution().seq() {
        println!("{}", x);
    }
    Ok(())
}
