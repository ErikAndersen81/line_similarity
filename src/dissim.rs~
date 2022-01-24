use geo::{prelude::HaversineDistance, Point};

pub fn similarity(trj_a: &[[f64; 3]], trj_b: &[[f64; 3]]) -> f64 {
    let timestamps: Vec<f64> = get_timestamps(trj_a, trj_b);
    let mut a: usize = 0;
    let mut b: usize = 0;
    let mut summands: Vec<f64> = vec![];
    assert!(timestamps.len() > 1);
    for k in 0..=timestamps.len() - 2 {
        let time = timestamps[k];
        let time_1 = timestamps[k + 1];
        if !(trj_a[a][2]..=trj_a[a + 1][2]).contains(&time) {
            a += 1;
        }
        if !(trj_b[b][2]..=trj_b[b + 1][2]).contains(&time) {
            b += 1;
        }
        let dist_1 = horizontal_euclidean(&trj_a[a], &trj_a[a + 1], &trj_b[b], &trj_b[b + 1], time);
        if !(trj_a[a][2]..=trj_a[a + 1][2]).contains(&time_1) {
            a += 1;
        }
        if !(trj_b[b][2]..=trj_b[b + 1][2]).contains(&time_1) {
            b += 1;
        }
        let dist_2 =
            horizontal_euclidean(&trj_a[a], &trj_a[a + 1], &trj_b[b], &trj_b[b + 1], time_1);
        summands.push((dist_1 + dist_2) * (time_1 - time));
    }
    summands.iter().sum::<f64>() / 2.0
}

fn get_timestamps(trj_a: &[[f64; 3]], trj_b: &[[f64; 3]]) -> Vec<f64> {
    let mut i: usize = 0;
    let mut j: usize = 0;
    let mut timestamps: Vec<f64> = vec![];
    while (i < trj_a.len()) && (j < trj_b.len()) {
        if trj_a[i][2] < trj_b[j][2] {
            timestamps.push(trj_a[i][2]);
            i += 1;
        } else {
            timestamps.push(trj_b[j][2]);
            j += 1;
        }
    }

    while i < trj_a.len() {
        timestamps.push(trj_a[i][2]);
        i += 1;
    }
    while j < trj_b.len() {
        timestamps.push(trj_b[j][2]);
        j += 1;
    }
    timestamps
}

fn horizontal_euclidean(p1: &[f64; 3], p2: &[f64; 3], q1: &[f64; 3], q2: &[f64; 3], t: f64) -> f64 {
    let t1: f64 = if p1[2] > q1[2] { p1[2] } else { q1[2] };
    let t2: f64 = if p2[2] < q2[2] { p2[2] } else { q2[2] };
    assert!(t1 <= t2);
    assert!((t1..=t2).contains(&t));
    let delta_t: f64 = (t - t1) / (t2 - t1);
    let q_x = q1[0] + (q2[0] - q1[0]) * delta_t;
    let p_x = p1[0] + (p2[0] - p1[0]) * delta_t;
    let q_y = q1[1] + (q2[1] - q1[1]) * delta_t;
    let p_y = p1[1] + (p2[1] - p1[1]) * delta_t;
    real_dist(&[q_x, q_y], &[p_x, p_y])
}

pub fn real_dist(a: &[f64; 2], b: &[f64; 2]) -> f64 {
    let p1 = Point::<f64>::from((a[1], a[0]));
    let p2 = Point::<f64>::from((b[1], b[0]));
    p1.haversine_distance(&p2)
}
