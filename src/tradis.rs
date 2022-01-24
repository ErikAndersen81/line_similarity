pub fn similarity(trj_a: &[[f64; 3]], trj_b: &[[f64; 3]]) -> f64 {
    // Note: Trajectories must be trimmed before calling
    assert!(trj_a.len() > 1);
    assert!(trj_b.len() > 1);
    let trj_ax: Vec<[f64; 2]> = trj_a.iter().map(|p| [p[0], p[2]]).collect();
    let trj_ay: Vec<[f64; 2]> = trj_a.iter().map(|p| [p[1], p[2]]).collect();
    let trj_bx: Vec<[f64; 2]> = trj_b.iter().map(|p| [p[0], p[2]]).collect();
    let trj_by: Vec<[f64; 2]> = trj_b.iter().map(|p| [p[1], p[2]]).collect();
    let area_x = area_of_polygon(trj_ax, trj_bx);
    let area_y = area_of_polygon(trj_ay, trj_by);
    area_x + area_y
}

fn area_of_polygon(trj_a: Vec<[f64; 2]>, trj_b: Vec<[f64; 2]>) -> f64 {
    let mut events = get_event_queue(trj_a, trj_b);
    let mut area = 0.0;
    assert!(events.len() >= 4);
    let initial_point_a = events.pop().unwrap().0;
    let initial_point_b = events.pop().unwrap().0;
    let mut event: ([f64; 2], bool) = events.pop().unwrap();
    let mut points: [[f64; 2]; 3] = [initial_point_a, initial_point_b, event.0];
    while !events.is_empty() {
        area += calc_area(&points[0], &points[1], &points[2]);
        if event.1 {
            let a = event.0;
            let b = events.pop().unwrap().0;
            event = events.pop().unwrap();
            points = [a, b, event.0];
        } else {
            event = events.pop().unwrap();
            points = [points[1], points[2], event.0];
        }
    }
    area += calc_area(&points[0], &points[1], &points[2]);
    area
}

fn get_event_queue(trj_a: Vec<[f64; 2]>, trj_b: Vec<[f64; 2]>) -> Vec<([f64; 2], bool)> {
    let mut events: Vec<([f64; 2], bool)> = Vec::new();
    let mut a_iter: std::iter::Peekable<std::slice::Iter<[f64; 2]>> = trj_a.iter().peekable();
    let mut b_iter: std::iter::Peekable<std::slice::Iter<[f64; 2]>> = trj_b.iter().peekable();
    let mut line_a: Line = Line {
        start: *a_iter.next().unwrap(),
        end: *a_iter.next().unwrap(),
    };
    let mut line_b: Line = Line {
        start: *b_iter.next().unwrap(),
        end: *b_iter.next().unwrap(),
    };

    while a_iter.peek().is_some() && b_iter.peek().is_some() {
        if let Some(intersection_point) = line_a.intersection(&line_b) {
            // Insert intersection event after two regular events
            if line_a.start[1] < line_b.start[1] {
                events.push((line_a.start, false));
                events.push((line_b.start, false));
            } else {
                events.push((line_b.start, false));
                events.push((line_a.start, false));
            }
            events.push((intersection_point, true));
            line_a.advance(a_iter.next().unwrap());
            line_b.advance(b_iter.next().unwrap());
        } else {
            // Insert regular event
            if line_a.start[1] < line_b.start[1] {
                events.push((line_a.start, false));
                line_a.advance(a_iter.next().unwrap());
            } else {
                events.push((line_b.start, false));
                line_b.advance(b_iter.next().unwrap());
            }
        }
    }

    fn fill_in_remaining(
        mut events: Vec<([f64; 2], bool)>,
        line_live: &mut Line,
        mut iter: std::iter::Peekable<std::slice::Iter<[f64; 2]>>,
        line_fixed: &Line,
    ) -> Vec<([f64; 2], bool)> {
        while iter.peek().is_some() {
            if let Some(intersection_point) = line_live.intersection(line_fixed) {
                events.push((line_live.start, false));
                events.push((intersection_point, true));
                line_live.advance(iter.next().unwrap());
            } else {
                events.push((line_live.start, false));
                line_live.advance(iter.next().unwrap());
            }
        }
        events.push((line_live.start, false));
        if let Some(intersection_point) = line_live.intersection(line_fixed) {
            events.push((intersection_point, true))
        }
        events
    }

    if a_iter.peek().is_some() {
        events.push((line_b.start, false));
        events = fill_in_remaining(events, &mut line_a, a_iter, &line_b);
    } else if b_iter.peek().is_some() {
        events.push((line_a.start, false));
        events = fill_in_remaining(events, &mut line_b, b_iter, &line_a);
    } else {
        if let Some(intersection_point) = line_a.intersection(&line_b) {
            events.push((intersection_point, true));
        }
        if line_a.start[1] < line_b.start[1] {
            events.push((line_a.start, false));
            events.push((line_b.start, false));
        } else {
            events.push((line_b.start, false));
            events.push((line_a.start, false));
        }
    }

    // insert the last points stored in line a and line b
    events.push((line_a.end, false));
    events.push((line_b.end, false));
    events
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
    fn get_events_queue_short_cross_test() {
        let trj_a = vec![[0., 0.], [1., 1.]];
        let trj_b = vec![[1., 0.], [0., 1.]];
        let events = get_event_queue(trj_a, trj_b);
        assert_eq!(events.len(), 5);
    }

    #[test]
    fn get_events_queue_long_cross_test() {
        let trj_a = vec![[0., 0.], [0., 6.]];
        let trj_b = vec![[1., 0.], [-1., 2.], [1., 4.], [-1., 6.]];
        let events = get_event_queue(trj_a, trj_b);
        assert_eq!(events.len(), 9);
        // Check placement of intersection points
        assert!(!events[0].1);
        assert!(!events[1].1);
        assert!(events[2].1);
        assert!(!events[3].1);
        assert!(events[4].1);
        assert!(!events[5].1);
        assert!(events[6].1);
        assert!(!events[7].1);
        assert!(!events[8].1);
    }

    #[test]
    fn get_events_queue_long_cross_test2() {
        let trj_a = vec![[0., 0.], [0., 2.], [0., 4.], [0., 6.]];
        let trj_b = vec![[1., 0.], [-1., 2.], [1., 4.], [-1., 6.]];
        let events = get_event_queue(trj_a, trj_b);
        assert_eq!(events.len(), 11);
    }
}
