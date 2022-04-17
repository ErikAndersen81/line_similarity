/// Determine the TRADIS distance for two 3D trajectories.
/// Both `trj_a` and `trj_b` must have time as the third dimension, i.e.
/// it must be monotone. The trajectories should also be trimmed s.t.
/// cover a common timespan, that is, `trj_a[0][2]==trj_b[0][2]` and
/// similarly for the last two points of each trajectory.
pub fn similarity(trj_a: &[[f64; 3]], trj_b: &[[f64; 3]]) -> f64 {
    // Note: Trajectories must be trimmed before calling
    assert!(trj_a.len() > 1);
    assert!(trj_b.len() > 1);
    let trj_ax: Vec<[f64; 2]> = trj_a.iter().map(|p| [p[0], p[2]]).collect();
    let trj_ay: Vec<[f64; 2]> = trj_a.iter().map(|p| [p[1], p[2]]).collect();
    let trj_bx: Vec<[f64; 2]> = trj_b.iter().map(|p| [p[0], p[2]]).collect();
    let trj_by: Vec<[f64; 2]> = trj_b.iter().map(|p| [p[1], p[2]]).collect();
    let area_x: f64 = get_polygons(&trj_ax, &trj_bx)
        .iter()
        .map(|polygon| triangulate_area(polygon.to_vec()))
        .sum();
    let area_y: f64 = get_polygons(&trj_ay, &trj_by)
        .iter()
        .map(|polygon| triangulate_area(polygon.to_vec()))
        .sum();
    area_x + area_y
}

/// Calculate the area of a t-monotone polygon
fn triangulate_area(mut polygon: Vec<[f64; 2]>) -> f64 {
    let mut area = 0.0;
    let mut last_pt = None;
    let mut left = polygon.remove(0);
    let mut right = polygon.pop().unwrap();
    while polygon.len() > 1 {
        let last = polygon.len() - 1;
        if polygon[0][1] < polygon[last][1] {
            last_pt = Some(left);
            left = polygon.remove(0);
        } else if polygon[last][1] < polygon[0][1] {
            last_pt = Some(right);
            right = polygon.pop().unwrap();
        } else if polygon[0][0] < polygon[last][0] {
            last_pt = Some(left);
            left = polygon.remove(0);
        } else if polygon[last][0] < polygon[0][0] {
            last_pt = Some(right);
            right = polygon.pop().unwrap();
        }
        if let Some(last_pt) = last_pt {
            area += calc_area(&left, &right, &last_pt);
        }
    }
    area
}

/// Returns a list of y-monotone polygons in the form of adjacency lists.
fn get_polygons(polyline_p: &[[f64; 2]], polyline_q: &[[f64; 2]]) -> Vec<Vec<[f64; 2]>> {
    let mut polygon_collection = vec![];
    let mut polygon = vec![];
    let mut p = 1;
    let mut q = 1;
    let mut line_p = Line {
        start: polyline_p[p - 1],
        end: polyline_p[p],
    };
    let mut line_q = Line {
        start: polyline_q[q - 1],
        end: polyline_q[q],
    };
    while (p < polyline_p.len()) || (q < polyline_q.len()) {
        if let Some(ints) = line_p.intersection(&line_q) {
            // Add final points to the polygon
            polygon.push(line_p.start);
            polygon.push(ints);
            polygon.insert(0, line_q.copy_start());
            // Adjust adjacency list to start with lowest t-value
            while polygon[0][1] > polygon[1][1] {
                let move_to_back = polygon.remove(0);
                polygon.push(move_to_back);
            }
            // Remove any consecutive duplicates
            polygon.dedup();
            polygon_collection.push(polygon);
            // Start building a new polygon
            polygon = vec![ints];
            line_p.start = ints;
            line_q.start = ints;
        } else {
            if line_p.is_lower_than(&line_q) {
                polygon.push(line_p.start);
                p += 1;
                if p < polyline_p.len() {
                    line_p.advance(&polyline_p[p]);
                }
            } else if line_q.is_lower_than(&line_p) {
                polygon.insert(0, line_q.start);
                q += 1;
                if q < polyline_q.len() {
                    line_q.advance(&polyline_q[q]);
                }
            } else {
                // both line segments ends at the same t-value, so we
                // can surely advance both. We could have done this
                // using non-exclusive if's but this is more readable.
                polygon.push(line_p.start);
                polygon.insert(0, line_q.start);
                p += 1;
                q += 1;
                if p < polyline_p.len() {
                    line_p.advance(&polyline_p[p]);
                }
                if q < polyline_q.len() {
                    line_q.advance(&polyline_q[q]);
                }
            }
        }
    }
    // Add final points to the polygon
    polygon.push(line_p.start);
    polygon.push(line_p.end);
    polygon.insert(0, line_q.copy_start());
    polygon.insert(0, line_q.copy_end());
    // Adjust adjacency list to start with lowest t-value
    while polygon[0][1] > polygon[1][1] {
        let move_to_back = polygon.remove(0);
        polygon.push(move_to_back);
    }
    // Remove any consecutive duplicates
    polygon.dedup();
    // Check for the odd case where trajectories end at an intersection point
    if polygon.len() > 1 {
        polygon_collection.push(polygon);
    }
    polygon_collection
}

fn calc_area(p: &[f64; 2], q: &[f64; 2], r: &[f64; 2]) -> f64 {
    0.5 * (p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1])).abs()
}

struct Line {
    start: [f64; 2],
    end: [f64; 2],
}

impl Line {
    fn intersection(&self, other: &Self) -> Option<[f64; 2]> {
        #[allow(clippy::float_cmp)]
        if (self.start[0] == other.start[0]) && (self.start[1] == other.start[1]) {
            // Special case where the lines start at the same point does not
            // count as an intersection.
            return None;
        }
        let (x1, y1, x2, y2) = (self.start[0], self.start[1], self.end[0], self.end[1]);
        let (x3, y3, x4, y4) = (other.start[0], other.start[1], other.end[0], other.end[1]);
        let denom: f64 = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        let t: f64 = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
        let u: f64 = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denom;
        if (0.0..=1.0).contains(&t) {
            return Some([x1 + t * (x2 - x1), y1 + t * (y2 - y1)]);
        } else if (0.0..=1.0).contains(&u) {
            return Some([x3 + u * (x4 - x3), y3 + u * (y4 - y3)]);
        };
        None
    }

    fn advance(&mut self, next: &[f64; 2]) {
        self.start = self.end;
        self.end = *next;
    }

    fn is_lower_than(&self, other: &Self) -> bool {
        self.start[1] < other.start[1]
    }

    fn copy_start(&self) -> [f64; 2] {
        self.start
    }
    fn copy_end(&self) -> [f64; 2] {
        self.end
    }
}

#[cfg(test)]
mod tradis_test {
    use super::*;

    #[test]
    fn intersect_multiple() {
        let polyline_a = vec![[0., 0.], [2., 2.], [0., 4.], [2., 6.]];
        let polyline_b = vec![[1., 0.], [1., 6.]];
        let expected_intersections = vec![[1., 1.], [1., 3.], [1., 5.]];
        let intersections = get_intersections(&polyline_a, &polyline_b);
        println!("intersections: {:?}", intersections);
        println!("expected: {:?}", expected_intersections);
        intersections
            .iter()
            .zip(expected_intersections.iter())
            .for_each(|(p, q)| assert!(equal(*p, *q)));
    }

    #[test]
    fn intersect_at_points_on_polyline() {
        let polyline_a = vec![[0., 0.], [1., 1.], [0., 2.]];
        let polyline_b = vec![[0., 0.], [-1., 1.], [0., 2.]];
        let expected_intersections = vec![[0., 2.]];
        let intersections = get_intersections(&polyline_a, &polyline_b);
        println!("intersections: {:?}", intersections);
        println!("expected: {:?}", expected_intersections);
        assert_eq!(intersections.len(), expected_intersections.len());
        intersections
            .iter()
            .zip(expected_intersections.iter())
            .for_each(|(p, q)| assert!(equal(*p, *q)));
    }

    #[test]
    fn cross_intersect() {
        let polyline_a = vec![[0., 0.], [2., 2.]];
        let polyline_b = vec![[2., 0.], [0., 2.]];
        let expected_intersections = vec![[1., 1.]];
        let intersections = get_intersections(&polyline_a, &polyline_b);
        println!("intersections: {:?}", intersections);
        println!("expected: {:?}", expected_intersections);
        assert_eq!(intersections.len(), expected_intersections.len());
        intersections
            .iter()
            .zip(expected_intersections.iter())
            .for_each(|(p, q)| assert!(equal(*p, *q)));
    }

    #[test]
    fn diamond_intersect() {
        let polyline_a = vec![[0., 0.], [1., 1.], [0., 2.], [1., 3.], [0., 4.]];
        let polyline_b = vec![[0., 0.], [-1., 1.], [0., 2.], [-1., 3.], [0., 4.]];
        let expected_intersections = vec![[0., 2.], [0., 4.]];
        let intersections = get_intersections(&polyline_a, &polyline_b);
        println!("intersections: {:?}", intersections);
        println!("expected: {:?}", expected_intersections);
        assert_eq!(intersections.len(), expected_intersections.len());
        intersections
            .iter()
            .zip(expected_intersections.iter())
            .for_each(|(p, q)| assert!(equal(*p, *q)));
    }

    #[test]
    fn cross_polygon() {
        let polyline_a = vec![[0., 0.], [2., 2.]];
        let polyline_b = vec![[2., 0.], [0., 2.]];
        let expected_poly_1 = vec![[2.0, 0.0], [0.0, 0.0], [1.0, 1.0]];
        let expected_poly_2 = vec![[1.0, 1.0], [2.0, 2.0], [0.0, 2.0]];
        let polygons = get_polygons(&polyline_a, &polyline_b);
        assert_eq!(polygons.len(), 2);
        println!("polygon_1: {:?}", polygons[0]);
        println!("expected: {:?}", expected_poly_1);
        println!("polygon_2: {:?}", polygons[1]);
        println!("expected: {:?}", expected_poly_2);
        assert_eq!(polygons.len(), 2);
        for (p, q) in expected_poly_1.iter().zip(polygons[0].clone()) {
            assert!(equal(*p, q));
        }
        for (p, q) in expected_poly_2.iter().zip(polygons[1].clone()) {
            assert!(equal(*p, q));
        }
    }

    #[test]
    fn diamonds_polygon() {
        let polyline_a = vec![[0., 0.], [1., 1.], [0., 2.], [1., 3.], [0., 4.]];
        let polyline_b = vec![[0., 0.], [-1., 1.], [0., 2.], [-1., 3.], [0., 4.]];
        let expected_poly_1 = vec![[0.0, 0.0], [1.0, 1.0], [0.0, 2.0], [-1.0, 1.0]];
        let expected_poly_2 = vec![[0.0, 2.0], [1.0, 3.0], [0.0, 4.0], [-1.0, 3.0]];
        let polygons = get_polygons(&polyline_a, &polyline_b);
        assert_eq!(polygons.len(), 2);
        println!("polygon_1: {:?}", polygons[0]);
        println!("expected: {:?}", expected_poly_1);
        println!("polygon_2: {:?}", polygons[1]);
        println!("expected: {:?}", expected_poly_2);
        assert_eq!(polygons.len(), 2);
        for (p, q) in expected_poly_1.iter().zip(polygons[0].clone()) {
            assert!(equal(*p, q));
        }
        for (p, q) in expected_poly_2.iter().zip(polygons[1].clone()) {
            assert!(equal(*p, q));
        }
    }

    #[test]
    fn diamonds_trangulate() {
        let diamond = vec![[0.0, 0.0], [1.0, 1.0], [0.0, 2.0], [-1.0, 1.0]];
        let area = triangulate_area(diamond);
        assert!((area - 2.0).abs() > 0.001);
    }

    #[test]
    fn jagged_trangulate() {
        let polygon = vec![
            [0.0, 0.0],
            [2.0, 1.0],
            [1.0, 2.0],
            [3.0, 4.0],
            [1.0, 5.0],
            [2.0, 6.0],
            [0.0, 7.0],
        ];
        let area = triangulate_area(polygon);
        assert!((area - 11.0).abs() > 0.001);
    }
}
