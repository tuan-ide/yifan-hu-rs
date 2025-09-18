use crate::vector::Vec2;

#[derive(Debug)]
struct QuadNode {
    center: Vec2,
    half_size: f64,
    mass: f64,
    center_of_mass: Vec2,
    children: [Option<usize>; 4],
    points: Vec<(usize, Vec2, f64)>,
    is_leaf: bool,
}

#[derive(Debug)]
pub struct BarnesHutTree {
    nodes: Vec<QuadNode>,
    root: usize,
    max_depth: usize,
    capacity: usize,
}

impl BarnesHutTree {
    pub fn new(points: &[(usize, Vec2, f64)], max_depth: usize, capacity: usize) -> Self {
        let mut min_x = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;

        for &(_, pos, _) in points {
            min_x = min_x.min(pos.x);
            max_x = max_x.max(pos.x);
            min_y = min_y.min(pos.y);
            max_y = max_y.max(pos.y);
        }

        if !min_x.is_finite() {
            min_x = -1.0;
            max_x = 1.0;
            min_y = -1.0;
            max_y = 1.0;
        }

        let span_x = (max_x - min_x).max(1e-6);
        let span_y = (max_y - min_y).max(1e-6);
        let half_size = 0.5 * span_x.max(span_y);
        let center = Vec2::new((min_x + max_x) * 0.5, (min_y + max_y) * 0.5);

        let mut tree = BarnesHutTree {
            nodes: Vec::new(),
            root: 0,
            max_depth: max_depth.max(1),
            capacity: capacity.max(1),
        };

        tree.root = tree.add_node(center, half_size);
        for &(idx, pos, weight) in points {
            tree.insert_point(tree.root, idx, pos, weight, 0);
        }
        tree.update_mass_properties(tree.root);
        tree
    }

    fn add_node(&mut self, center: Vec2, half_size: f64) -> usize {
        let node = QuadNode {
            center,
            half_size,
            mass: 0.0,
            center_of_mass: Vec2::zero(),
            children: [None, None, None, None],
            points: Vec::new(),
            is_leaf: true,
        };
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    fn insert_point(&mut self, node_idx: usize, idx: usize, pos: Vec2, weight: f64, depth: usize) {
        if self.nodes[node_idx].is_leaf {
            if self.nodes[node_idx].points.len() < self.capacity || depth >= self.max_depth {
                self.nodes[node_idx].points.push((idx, pos, weight));
                return;
            }
            self.subdivide(node_idx, depth);
        }

        // Determine quadrant and insert into child.
        let quadrant = self.select_quadrant(node_idx, pos);
        let child_idx = self.nodes[node_idx].children[quadrant].unwrap();
        self.insert_point(child_idx, idx, pos, weight, depth + 1);
    }

    fn subdivide(&mut self, node_idx: usize, depth: usize) {
        let (center, half_size);
        {
            let node = &mut self.nodes[node_idx];
            node.is_leaf = false;
            center = node.center;
            half_size = node.half_size * 0.5;
        }

        let offsets = [
            Vec2::new(-0.5, -0.5),
            Vec2::new(0.5, -0.5),
            Vec2::new(-0.5, 0.5),
            Vec2::new(0.5, 0.5),
        ];

        for (i, offset) in offsets.iter().enumerate() {
            let child_center = Vec2::new(
                center.x + offset.x * 2.0 * half_size,
                center.y + offset.y * 2.0 * half_size,
            );
            let child_idx = self.add_node(child_center, half_size);
            self.nodes[node_idx].children[i] = Some(child_idx);
        }

        let points = std::mem::take(&mut self.nodes[node_idx].points);
        for (idx, pos, weight) in points {
            let quadrant = self.select_quadrant(node_idx, pos);
            let child_idx = self.nodes[node_idx].children[quadrant].unwrap();
            self.insert_point(child_idx, idx, pos, weight, depth + 1);
        }
    }

    fn select_quadrant(&self, node_idx: usize, pos: Vec2) -> usize {
        let node = &self.nodes[node_idx];
        let east = if pos.x >= node.center.x { 1 } else { 0 };
        let north = if pos.y >= node.center.y { 2 } else { 0 };
        east + north
    }

    fn update_mass_properties(&mut self, node_idx: usize) -> (f64, Vec2) {
        let (mass, com) = if self.nodes[node_idx].is_leaf {
            let mut mass = 0.0;
            let mut weighted = Vec2::zero();
            for &(_, pos, weight) in &self.nodes[node_idx].points {
                mass += weight;
                weighted += pos * weight;
            }
            if mass > 0.0 {
                (mass, weighted / mass)
            } else {
                (0.0, self.nodes[node_idx].center)
            }
        } else {
            let mut mass = 0.0;
            let mut weighted = Vec2::zero();
            let children: Vec<usize> = self.nodes[node_idx]
                .children
                .iter()
                .flatten()
                .copied()
                .collect();
            for child_idx in children {
                let (child_mass, child_com) = self.update_mass_properties(child_idx);
                if child_mass > 0.0 {
                    mass += child_mass;
                    weighted += child_com * child_mass;
                }
            }
            if mass > 0.0 {
                (mass, weighted / mass)
            } else {
                (0.0, self.nodes[node_idx].center)
            }
        };
        self.nodes[node_idx].mass = mass;
        self.nodes[node_idx].center_of_mass = com;
        (mass, com)
    }

    pub fn repulsive_force(
        &self,
        index: usize,
        position: Vec2,
        theta: f64,
        cutoff: Option<f64>,
        coeff: f64,
        exponent: f64,
    ) -> Vec2 {
        if self.nodes.is_empty() {
            return Vec2::zero();
        }
        let mut force = Vec2::zero();
        let mut stack = Vec::new();
        stack.push(self.root);

        while let Some(node_idx) = stack.pop() {
            let node = &self.nodes[node_idx];
            if node.mass == 0.0 {
                continue;
            }
            if node.is_leaf {
                for &(other_idx, other_pos, other_weight) in &node.points {
                    if other_idx == index {
                        continue;
                    }
                    force +=
                        repulsive_pair(position, other_pos, other_weight, cutoff, coeff, exponent);
                }
            } else {
                let diff = position - node.center_of_mass;
                let dist_sq = diff.length_squared();
                if dist_sq == 0.0 {
                    // If the vertex lies exactly at the center of mass, descend to children to avoid division by zero.
                    for child in node.children.iter().flatten() {
                        stack.push(*child);
                    }
                    continue;
                }
                let dist = dist_sq.sqrt();
                if cutoff.is_some() && dist > cutoff.unwrap() {
                    continue;
                }
                let width = node.half_size * 2.0;
                if width / dist < theta {
                    let scaled = diff * (coeff * node.mass / dist.powf(exponent + 1.0));
                    force += scaled;
                } else {
                    for child in node.children.iter().flatten() {
                        stack.push(*child);
                    }
                }
            }
        }

        force
    }
}

fn repulsive_pair(
    from: Vec2,
    other: Vec2,
    other_weight: f64,
    cutoff: Option<f64>,
    coeff: f64,
    exponent: f64,
) -> Vec2 {
    let delta = from - other;
    let dist_sq = delta.length_squared();
    if dist_sq == 0.0 {
        return Vec2::zero();
    }
    let dist = dist_sq.sqrt();
    if let Some(cutoff) = cutoff {
        if dist > cutoff {
            return Vec2::zero();
        }
    }
    let scale = coeff * other_weight / dist.powf(exponent + 1.0);
    delta * scale
}
