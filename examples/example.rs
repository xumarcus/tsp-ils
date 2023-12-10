use std::{
    collections::HashMap,
    ffi::OsStr,
    fs::{self, File},
    io,
    path::Path,
    time::Duration,
};

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use statrs::statistics::{Data, Distribution, Max, Min};
use tsp_ils::{Solver, Tour, TSP};

fn main() -> io::Result<()> {
    let optimal = HashMap::from([
        (OsStr::new("d657.tsp"), 48912),
        (OsStr::new("dsj1000.tsp"), 18659688),
        (OsStr::new("lu980.tsp"), 11340),
        (OsStr::new("uy734.tsp"), 79114),
        (OsStr::new("zi929.tsp"), 95345),
    ]);

    for entry in fs::read_dir(Path::new("examples/data"))? {
        let path = entry?.path();
        let file = File::open(path.clone())?;
        let tsp = TSP::parse_tsp_file(file)?;

        let filename = path.file_name().expect(".. not allowed");
        let opt = *optimal.get(filename).expect("no optimal recorded") as f64;
        let naive = Tour::new(&tsp).cost() as f64;
        let results: Vec<f64> = (0..12)
            .into_par_iter()
            .map(|seed| {
                let mut s = Solver::new(&tsp, Duration::from_secs(2), seed);
                s.run();
                let val = s.solution().cost() as f64;
                0.02_f64.powf((val - opt) / (naive - val))
            })
            .collect();
        let data = Data::new(results);

        println!(
            "{:>16?}: min={:.03} max={:.03} mean={:.03} stdev={:.03}",
            filename,
            data.min(),
            data.max(),
            data.mean().unwrap(),
            data.std_dev().unwrap()
        );
    }
    Ok(())
}
