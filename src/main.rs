use std::{io, time::Duration};

use tsp_ils::{Solver, TSP};

fn main() -> io::Result<()> {
    let t = TSP::parse_kattis()?;
    let mut s = Solver::new(&t, Duration::from_secs(2), 4234);
    s.run();
    for x in s.solution().seq() {
        println!("{}", x);
    }
    Ok(())
}
