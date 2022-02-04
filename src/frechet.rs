use crate::euclidean_distance;

pub fn similarity(trj_a: &[[f64; 3]], trj_b: &[[f64; 3]]) -> f64 {
    // Note: Trajectories must be trimmed before calling
    assert!(trj_a.len() > 1);
    assert!(trj_b.len() > 1);
    let trj_a: Vec<[f64; 2]> = trj_a.iter().map(|[x, y, _]| [*x, *y]).collect();
    let trj_b: Vec<[f64; 2]> = trj_b.iter().map(|[x, y, _]| [*x, *y]).collect();
    frechet_distance(&trj_a, &trj_b)
}

fn frechet_distance(ls_a: &[[f64; 2]], ls_b: &[[f64; 2]]) -> f64 {
    let mut data = Data {
        cache: vec![vec![f64::NAN; ls_b.len()]; ls_a.len()],
        ls_a,
        ls_b,
    };
    data.compute(ls_a.len() - 1, ls_b.len() - 1)
}

struct Data<'a> {
    cache: Vec<Vec<f64>>,
    ls_a: &'a [[f64; 2]],
    ls_b: &'a [[f64; 2]],
}

impl<'a> Data<'a> {
    fn compute(&mut self, i: usize, j: usize) -> f64 {
        if self.cache[i][j].is_nan() {
            let dist = euclidean_distance(&self.ls_a[i], &self.ls_b[j]);
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
