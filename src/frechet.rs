use geo::coords_iter::CoordsIter;
use geo::prelude::HaversineDistance;
use geo::LineString;
use geo::Point;

pub fn similarity(trj_a: &[[f64; 3]], trj_b: &[[f64; 3]]) -> f64 {
    // Note: Trajectories must be trimmed before calling
    assert!(trj_a.len() > 1);
    assert!(trj_b.len() > 1);
    let trj_a = trj_a
        .iter()
        .map(|p| Point::<f64>::from((p[0], p[1])))
        .collect::<Vec<Point<f64>>>();
    let trj_b: Vec<Point<f64>> = trj_b
        .iter()
        .rev()
        .map(|p| Point::<f64>::from((p[0], p[1])))
        .collect();
    let trj_a = LineString::from(trj_a);
    let trj_b = LineString::from(trj_b);
    frechet_distance(&trj_a, &trj_b)
}

fn frechet_distance(ls_a: &LineString<f64>, ls_b: &LineString<f64>) -> f64 {
    if ls_a.coords_count() != 0 && ls_b.coords_count() != 0 {
        let mut data = Data {
            cache: vec![vec![f64::NAN; ls_b.coords_count()]; ls_a.coords_count()],
            ls_a,
            ls_b,
        };
        data.compute(ls_a.coords_count() - 1, ls_b.coords_count() - 1)
    } else {
        0.0
    }
}

struct Data<'a> {
    cache: Vec<Vec<f64>>,
    ls_a: &'a LineString<f64>,
    ls_b: &'a LineString<f64>,
}

impl<'a> Data<'a> {
    fn compute(&mut self, i: usize, j: usize) -> f64 {
        if self.cache[i][j].is_nan() {
            let dist = Point::from(self.ls_a[i]).haversine_distance(&Point::from(self.ls_b[j]));
            self.cache[i][j] = match (i, j) {
                (0, 0) => dist,
                (_, 0) => self.compute(i - 1, 0).max(dist),
                (0, _) => self.compute(0, j - 1).max(dist),
                (_, _) => ((self.compute(i - 1, j).min(self.compute(i - 1, j - 1)))
                    .min(self.compute(i, j - 1)))
                .max(dist),
            };
        }
        self.cache[i][j]
    }
}
