use crate::euclidean_distance;

pub fn similarity(trj_a: &[[f64; 3]], trj_b: &[[f64; 3]]) -> f64 {
    let mut infima_a = vec![f64::NAN; trj_a.len()];
    let mut infima_b = vec![f64::NAN; trj_b.len()];
    for i in 0..trj_a.len() {
        for j in 0..trj_b.len() {
            let a = [trj_a[i][0], trj_a[i][1]];
            let b = [trj_b[i][0], trj_b[i][1]];
            let dist = euclidean_distance(&a, &b);
            infima_a[i] = infima_a[i].min(dist);
            infima_b[j] = infima_b[j].min(dist);
        }
    }
    let sup_inf_a: f64 = infima_a.iter().fold(f64::NAN, |acc, val| acc.max(*val));
    let sup_inf_b: f64 = infima_b.iter().fold(f64::NAN, |acc, val| acc.max(*val));
    sup_inf_a.max(sup_inf_b)
}
