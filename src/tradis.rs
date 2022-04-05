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

/// Calculate the area of a y-monotone polygon
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
fn get_polygons(polyline_a: &[[f64; 2]], polyline_b: &[[f64; 2]]) -> Vec<Vec<[f64; 2]>> {
    let intersections = get_intersections(polyline_a, polyline_b);
    // Remove any points from the polylines that are also intersection points.
    let polyline_a: Vec<[f64; 2]> = polyline_a
        .iter()
        .filter(|&p| intersections.iter().all(|&q| !equal(*p, q)))
        .copied()
        .collect();
    let polyline_b: Vec<[f64; 2]> = polyline_b
        .iter()
        .filter(|&p| intersections.iter().all(|&q| !equal(*p, q)))
        .copied()
        .collect();

    let mut polygon_collection: Vec<Vec<[f64; 2]>> = vec![];
    #[derive(Copy, Clone)]
    enum Case {
        /// Single intersection point at top of polygon
        Bottom,
        /// Intersection points at top and bottom of polygon
        Middle,
        /// Single intersection point at bottom of polygon
        Top,
    }
    let mut last_case = None;
    for (idx, intersection) in intersections.iter().enumerate() {
        if idx + 1 == intersections.len() && intersection[1] > polyline_a[polyline_a.len() - 1][1] {
            // Special case where the polygon ends with a Case::Middle
            last_case = Some(Case::Top);
            continue;
        }
        let case: Case = if idx == 0 && intersection[1] > polyline_a[0][1] {
            Case::Bottom
        } else if idx + 1 == intersections.len()
            && intersection[1] < polyline_a[polyline_a.len() - 1][1]
        {
            Case::Top
        } else {
            Case::Middle
        };
        last_case = Some(case);
        let interval = |case| match case {
            Case::Bottom => 0.0..intersection[1],
            Case::Top => intersection[1]..f64::INFINITY,
            Case::Middle => intersections[idx][1]..intersections[idx + 1][1],
        };
        let polygon = match case {
            Case::Top => {
                let mut polygon = vec![*intersection];
                let points_a = polyline_a.iter().filter(|p| interval(case).contains(&p[1]));
                polygon.extend(points_a);
                let points_b = polyline_b
                    .iter()
                    .filter(|p| interval(case).contains(&p[1]))
                    .rev();
                polygon.extend(points_b);
                polygon
            }
            Case::Middle => {
                let mut polygon = vec![*intersection];
                let points_a = polyline_a.iter().filter(|p| interval(case).contains(&p[1]));
                polygon.extend(points_a);
                polygon.push(intersections[idx + 1]);
                let points_b = polyline_b
                    .iter()
                    .filter(|p| interval(case).contains(&p[1]))
                    .rev();
                polygon.extend(points_b);
                polygon
            }
            Case::Bottom => {
                let mut polygon = polyline_a
                    .iter()
                    .filter(|p| interval(case).contains(&p[1]))
                    .copied()
                    .collect::<Vec<[f64; 2]>>();
                polygon.push(*intersection);
                let points_b = polyline_b
                    .iter()
                    .filter(|p| interval(case).contains(&p[1]))
                    .rev();
                polygon.extend(points_b);
                polygon
            }
        };
        polygon_collection.push(polygon);
    }
    if let Some(case) = last_case {
        match case {
            Case::Middle | Case::Bottom => {
                // Construct the final polygon
                let intersection = intersections[intersections.len() - 1];
                let mut polygon = vec![intersection];
                let points_a = polyline_a.iter().filter(|p| intersection[1] < p[1]);
                polygon.extend(points_a);
                let points_b = polyline_b.iter().filter(|p| intersection[1] < p[1]).rev();
                polygon.extend(points_b);
                polygon_collection.push(polygon);
            }
            _ => {}
        }
    }
    polygon_collection
}

fn get_intersections(polyline_a: &[[f64; 2]], polyline_b: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut intersections: Vec<[f64; 2]> = vec![];
    let mut i: usize = 1;
    let mut j: usize = 1;
    let mut line_a: Line = Line {
        start: polyline_a[i - 1],
        end: polyline_a[i],
    };
    let mut line_b: Line = Line {
        start: polyline_b[j - 1],
        end: polyline_b[j],
    };
    while (i < polyline_a.len()) || (j < polyline_b.len()) {
        if let Some(intersection_point) = line_a.intersection(&line_b) {
            if intersections.is_empty()
                || !equal(intersection_point, intersections[intersections.len() - 1])
            {
                // Handle the special case where intersection points are coincidental with polyline points
                intersections.push(intersection_point);
            }
        }
        // Advance line with lowest y-value on end-coordinate
        if line_a.end[1] < line_b.end[1] {
            i += 1;
            if i < polyline_a.len() {
                line_a.advance(&polyline_a[i]);
            }
        } else if line_b.end[1] < line_a.end[1] {
            j += 1;
            if j < polyline_b.len() {
                line_b.advance(&polyline_b[j]);
            }
        } else {
            // Advance line with highest y-value on start-coordinate
            j += 1;
            i += 1;
            if j < polyline_b.len() {
                line_b.advance(&polyline_b[j]);
            }
            if i < polyline_a.len() {
                line_a.advance(&polyline_a[i]);
            }
        }
    }
    intersections
}

/// Determine if two 2D-points are similar within epsilon 0.00000001
fn equal(p: [f64; 2], q: [f64; 2]) -> bool {
    (p[0] - q[0]).abs() + (p[1] - q[1]).abs() < 0.00000001
}

fn calc_area(p: &[f64; 2], q: &[f64; 2], r: &[f64; 2]) -> f64 {
    0.5 * (p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1])).abs()
}

struct Line {
    start: [f64; 2],
    end: [f64; 2],
}

impl Line {
    fn intersection(&self, other: &Line) -> Option<[f64; 2]> {
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
        let expected_intersections = vec![[0., 0.], [0., 2.]];
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
        let expected_intersections = vec![[0., 0.], [0., 2.], [0., 4.]];
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
        let expected_poly_1 = vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]];
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
