use std::{
    collections::HashMap,
    ffi::OsStr,
    fs::{self, File},
    io::{self, Error, ErrorKind},
    path::Path,
    time::Duration,
};

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use statrs::statistics::{Data, Distribution, Max, Min};
use tsp_ils::{Solver, Tour, TSP};

fn summary<D: AsMut<[f64]> + AsRef<[f64]>>(data: Data<D>) -> Option<String> {
    Some(format!(
        "min={:.03} max={:.03} mean={:.03} stdev={:.03}",
        data.min(),
        data.max(),
        data.mean()?,
        data.std_dev()?
    ))
}

fn main() -> io::Result<()> {
    let optimal = HashMap::from([
        ("d657.tsp", 48912),
        ("dsj1000.tsp", 18659688),
        ("lu980.tsp", 11340),
        ("uy734.tsp", 79114),
        ("zi929.tsp", 95345),
    ]);

    let ms = fs::read_dir(Path::new("examples/data"))?
        .map(|entry| {
            let path = entry?.path();
            let file = File::open(path.clone())?;
            let tsp = TSP::parse_tsp_file(file)?;

            let filename = path.file_name().and_then(OsStr::to_str).expect("valid str");
            let opt = *optimal.get(filename).expect("no optimal recorded") as f64;
            let naive = Tour::new(&tsp).cost() as f64;
            let (cv, vv): (Vec<_>, Vec<_>) = (0..12)
                .into_par_iter()
                .map(|seed| {
                    let mut s = Solver::new(&tsp, Duration::from_secs(2), seed);
                    let cnt = s.run();
                    let val = s.solution().cost() as f64;
                    (
                        1e-3 * cnt as f64,
                        0.02_f64.powf((val - opt) / (naive - val)),
                    )
                })
                .unzip();
            let cd = Data::new(cv);
            let vd = Data::new(vv);
            let (vm, cs, vs) = (|| {
                let vm = vd.mean()?;
                let cs = summary(cd)?;
                let vs = summary(vd)?;
                Some((vm, cs, vs))
            })()
            .ok_or(Error::new(ErrorKind::Other, "Empty data"))?;
            println!("{:>16}: [cnt/1k] {} [val] {}", filename, cs, vs);
            Ok(vm)
        })
        .collect::<io::Result<Vec<f64>>>()?;
    let md = Data::new(ms);
    println!("suite: {}", summary(md).unwrap());
    Ok(())
}
