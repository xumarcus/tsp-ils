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
        ("d657.tsp", 48912),
        ("dsj1000.tsp", 18659688),
        ("lu980.tsp", 11340),
        ("uy734.tsp", 79114),
        ("zi929.tsp", 95345),
    ]);

    let mut means: Vec<f64> = Vec::new();
    for entry in fs::read_dir(Path::new("examples/data"))? {
        let path = entry?.path();
        let file = File::open(path.clone())?;
        let tsp = TSP::parse_tsp_file(file)?;

        let filename = path.file_name().and_then(OsStr::to_str).expect("valid str");
        let opt = *optimal.get(filename).expect("no optimal recorded") as f64;
        let naive = Tour::new(&tsp).cost() as f64;
        let (cnts, scores): (Vec<_>, Vec<_>) = (0..12)
            .into_par_iter()
            .map(|seed| {
                let mut s = Solver::new(&tsp, Duration::from_secs(2), seed);
                let cnt = s.run();
                let val = s.solution().cost() as f64;
                (cnt as f64, 0.02_f64.powf((val - opt) / (naive - val)))
            })
            .unzip();
        let cnt_data = Data::new(cnts);
        let score_data = Data::new(scores);
        means.push(score_data.mean().unwrap());

        println!(
            "{:>16}: [cnt] min={:>6} max={:>6} mean={:>6.0} stdev={:>6.0} [val] min={:.03} max={:.03} mean={:.03} stdev={:.03}",
            filename,
            cnt_data.min(),
            cnt_data.max(),
            cnt_data.mean().unwrap(),
            cnt_data.std_dev().unwrap(),
            score_data.min(),
            score_data.max(),
            score_data.mean().unwrap(),
            score_data.std_dev().unwrap()
        );
    }
    let mean_data = Data::new(means);
    println!(
        "mean: {:.03} stdev={:.03}",
        mean_data.mean().unwrap(),
        mean_data.std_dev().unwrap()
    );
    Ok(())
}
