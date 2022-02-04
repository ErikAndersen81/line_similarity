pub mod dissim;
pub mod frechet;
pub mod hausdorff;
pub mod tradis;

pub fn euclidean_distance(a: &[f64; 2], b: &[f64; 2]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}
