use std::io;

use tsp_ils::{Solver, TSP};

fn main() -> io::Result<()> {
    /*
    let t = File::open(Path::new("tests/sample.tsp"))
        .and_then(TSP::parse_tsp_file)?;
     */
    let t = TSP::parse_kattis()?;
    let mut s = Solver::new(&t);
    s.run();
    for x in s.solution().seq() {
        println!("{}", x);
    }
    Ok(())
}
