use std::{path::Path, io};

use tsp_ils::{TSP, Solver};

fn main() -> io::Result<()> {
    let t = TSP::from_path(Path::new("tests/berlin52.tsp"))?;
    let mut s = Solver::new(&t);
    // s.run();
    println!("{}", s.solution().cost());
    Ok(())
}
