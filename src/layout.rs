use crate::barnes_hut::BarnesHutTree;
use crate::coarsen::{CoarsenResult, CoarseningStrategy, coarsen};
use crate::graph::Graph;
use crate::vector::Vec2;
use rand::prelude::*;

#[derive(Debug, Clone)]
pub struct LayoutSettings {
    pub tolerance: f64,
    pub max_iterations: usize,
    pub initial_step: f64,
    pub min_step: f64,
    pub max_step: f64,
    pub cooling_factor: f64,
    pub adaptive_decay: f64,
    pub adaptive_progress_limit: usize,
    pub repulsive_exponent: f64,
    pub theta: f64,
    pub c_strength: f64,
    pub repulsive_cutoff: Option<f64>,
    pub tree_capacity: usize,
    pub tree_depth: usize,
    pub coarsening_strategy: CoarseningStrategy,
    pub coarsening_ratio_threshold: f64,
    pub min_coarse_size: usize,
    pub random_seed: Option<u64>,
    pub jitter_fraction: f64,
}

impl Default for LayoutSettings {
    fn default() -> Self {
        Self {
            tolerance: 1e-3,
            max_iterations: 500,
            initial_step: 1.0,
            min_step: 1e-4,
            max_step: 10.0,
            cooling_factor: 0.9,
            adaptive_decay: 0.9,
            adaptive_progress_limit: 5,
            repulsive_exponent: 1.0,
            theta: 1.2,
            c_strength: 0.2,
            repulsive_cutoff: None,
            tree_capacity: 1,
            tree_depth: 10,
            coarsening_strategy: CoarseningStrategy::Hybrid,
            coarsening_ratio_threshold: 0.75,
            min_coarse_size: 2,
            random_seed: None,
            jitter_fraction: 0.01,
        }
    }
}

#[derive(Debug)]
pub struct LayoutResult {
    pub positions: Vec<Vec2>,
    pub iterations: usize,
}

enum StepScheme {
    Adaptive,
    Simple,
}

struct ForceParams<'a> {
    graph: &'a Graph,
    positions: &'a mut [Vec2],
    step_scheme: StepScheme,
    settings: &'a LayoutSettings,
}

pub fn multilevel_layout(graph: &Graph, settings: &LayoutSettings) -> LayoutResult {
    let seed = settings.random_seed.unwrap_or(42);
    let mut rng = StdRng::seed_from_u64(seed);
    let mut positions = vec![Vec2::zero(); graph.node_count()];
    let iterations = multilevel_layout_recursive(graph, settings, &mut rng, &mut positions);
    LayoutResult {
        positions,
        iterations,
    }
}

fn multilevel_layout_recursive(
    graph: &Graph,
    settings: &LayoutSettings,
    rng: &mut StdRng,
    positions: &mut [Vec2],
) -> usize {
    if graph.node_count() <= settings.min_coarse_size {
        randomize_positions(positions, rng);
        return force_directed(ForceParams {
            graph,
            positions,
            step_scheme: StepScheme::Adaptive,
            settings,
        });
    }

    match coarsen(
        graph,
        settings.coarsening_strategy,
        settings.coarsening_ratio_threshold,
        rng,
    ) {
        Some(CoarsenResult {
            coarse,
            prolongation,
            ..
        }) if coarse.node_count() >= settings.min_coarse_size => {
            let mut coarse_positions = vec![Vec2::zero(); coarse.node_count()];
            let coarse_iterations =
                multilevel_layout_recursive(&coarse, settings, rng, &mut coarse_positions);

            prolongate(&prolongation, &coarse_positions, positions);
            jitter_overlaps(positions, rng, settings.jitter_fraction);

            let fine_diameter = graph.approximate_diameter().max(1) as f64;
            let coarse_diameter = coarse.approximate_diameter().max(1) as f64;
            if coarse_diameter > 0.0 {
                let scale = fine_diameter / coarse_diameter;
                for pos in positions.iter_mut() {
                    *pos = *pos * scale;
                }
            }

            let refine_iterations = force_directed(ForceParams {
                graph,
                positions,
                step_scheme: StepScheme::Simple,
                settings,
            });
            coarse_iterations + refine_iterations
        }
        _ => {
            randomize_positions(positions, rng);
            force_directed(ForceParams {
                graph,
                positions,
                step_scheme: StepScheme::Adaptive,
                settings,
            })
        }
    }
}

fn randomize_positions(positions: &mut [Vec2], rng: &mut StdRng) {
    let count = positions.len() as f64;
    let radius = (count.sqrt()).max(1.0);
    for pos in positions.iter_mut() {
        let x = rng.gen_range(-radius..radius);
        let y = rng.gen_range(-radius..radius);
        *pos = Vec2::new(x, y);
    }
}

fn jitter_overlaps(positions: &mut [Vec2], rng: &mut StdRng, fraction: f64) {
    if positions.is_empty() {
        return;
    }
    let mut lookup = std::collections::HashMap::<(i64, i64), Vec<usize>>::new();
    let scale = 1.0 / fraction.max(1e-6);
    for (idx, pos) in positions.iter().enumerate() {
        let key = (
            (pos.x * scale).round() as i64,
            (pos.y * scale).round() as i64,
        );
        lookup.entry(key).or_default().push(idx);
    }
    for group in lookup.values() {
        if group.len() <= 1 {
            continue;
        }
        for &idx in group {
            let jitter_x = rng.gen_range(-fraction..fraction);
            let jitter_y = rng.gen_range(-fraction..fraction);
            positions[idx].x += jitter_x;
            positions[idx].y += jitter_y;
        }
    }
}

pub(crate) fn update_step_adaptive(
    step: f64,
    progress: usize,
    energy: f64,
    energy_prev: f64,
    settings: &LayoutSettings,
) -> (f64, usize, f64) {
    let mut step = step;
    let mut progress = progress;
    if energy < energy_prev {
        progress += 1;
        if progress >= settings.adaptive_progress_limit {
            let denom = settings.adaptive_decay.max(1e-6);
            step = (step / denom).min(settings.max_step);
            progress = 0;
        }
    } else {
        step = (step * settings.adaptive_decay).max(settings.min_step);
        progress = 0;
    }
    (step, progress, energy)
}

fn prolongate(
    prolongation: &[Vec<(usize, f64)>],
    coarse_positions: &[Vec2],
    fine_positions: &mut [Vec2],
) {
    for (fine_idx, mapping) in prolongation.iter().enumerate() {
        if mapping.is_empty() {
            fine_positions[fine_idx] = Vec2::zero();
            continue;
        }
        let mut pos = Vec2::zero();
        for &(coarse_idx, weight) in mapping {
            pos += coarse_positions[coarse_idx] * weight;
        }
        fine_positions[fine_idx] = pos;
    }
}

fn force_directed(params: ForceParams<'_>) -> usize {
    let ForceParams {
        graph,
        positions,
        step_scheme,
        settings,
        ..
    } = params;

    let mut step = settings
        .initial_step
        .clamp(settings.min_step, settings.max_step);
    let mut energy_prev = f64::INFINITY;
    let mut progress = 0usize;
    let mut iterations = 0usize;

    let mut snapshot = positions.to_vec();

    while iterations < settings.max_iterations {
        snapshot.clone_from_slice(positions);
        let points: Vec<(usize, Vec2, f64)> = snapshot
            .iter()
            .enumerate()
            .map(|(i, pos)| (i, *pos, graph.node_weight(i)))
            .collect();
        let tree = BarnesHutTree::new(&points, settings.tree_depth, settings.tree_capacity);

        let avg_edge = average_edge_length(graph, &snapshot).max(1e-6);
        let coeff = settings.c_strength * avg_edge.powf(settings.repulsive_exponent + 1.0);

        let mut energy = 0.0;
        let mut displacement = 0.0;

        for i in 0..graph.node_count() {
            let mut force = Vec2::zero();
            let origin = snapshot[i];

            // Attractive forces.
            for neighbor in graph.neighbors(i).iter() {
                let target = neighbor.target;
                let weight = neighbor.weight;
                let delta = origin - positions[target];
                let dist = delta.length().max(1e-9);
                let attractive = delta * (-weight * dist / avg_edge);
                force += attractive;
            }

            // Repulsive forces via Barnes-Hut.
            let repulsive = tree.repulsive_force(
                i,
                positions[i],
                settings.theta,
                settings.repulsive_cutoff,
                coeff,
                settings.repulsive_exponent,
            );
            force += repulsive;

            let magnitude = force.length();
            if magnitude > 0.0 {
                let direction = force / magnitude;
                let limited_step = step.min(settings.max_step);
                positions[i] += direction * limited_step;
                displacement += (positions[i] - origin).length();
            }
            energy += magnitude * magnitude;
        }

        iterations += 1;
        let avg_displacement = displacement / graph.node_count() as f64;
        if avg_displacement < settings.tolerance {
            break;
        }

        match step_scheme {
            StepScheme::Adaptive => {
                let (new_step, new_progress, new_energy_prev) =
                    update_step_adaptive(step, progress, energy, energy_prev, settings);
                step = new_step;
                progress = new_progress;
                energy_prev = new_energy_prev;
            }
            StepScheme::Simple => {
                step = (step * settings.cooling_factor).clamp(settings.min_step, settings.max_step);
                energy_prev = energy;
            }
        }
    }

    iterations
}

fn average_edge_length(graph: &Graph, positions: &[Vec2]) -> f64 {
    let mut total = 0.0;
    let mut count = 0.0;
    for (u, node) in graph.nodes.iter().enumerate() {
        for neighbor in &node.neighbors {
            if neighbor.target < u {
                continue;
            }
            let delta = positions[u] - positions[neighbor.target];
            total += delta.length();
            count += 1.0;
        }
    }
    if count > 0.0 { total / count } else { 1.0 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;

    #[test]
    fn small_graph_layout_is_finite() {
        let mut graph = Graph::new(4);
        graph.add_edge(0, 1, 1.0);
        graph.add_edge(1, 2, 1.0);
        graph.add_edge(2, 3, 1.0);
        let mut settings = LayoutSettings::default();
        settings.max_iterations = 50;
        settings.tolerance = 1e-2;
        settings.random_seed = Some(123);
        let result = multilevel_layout(&graph, &settings);
        assert_eq!(result.positions.len(), 4);
        for pos in result.positions {
            assert!(pos.x.is_finite() && pos.y.is_finite());
        }
    }

    #[test]
    fn adaptive_step_matches_pdf_progress_increase() {
        let settings = LayoutSettings::default();
        let (step, progress, energy_prev) = update_step_adaptive(1.0, 4, 0.5, 1.0, &settings);
        let expected = 1.0 / settings.adaptive_decay;
        assert!((step - expected).abs() < 1e-9);
        assert_eq!(progress, 0);
        assert!((energy_prev - 0.5).abs() < 1e-12);
    }

    #[test]
    fn adaptive_step_matches_pdf_energy_increase() {
        let settings = LayoutSettings::default();
        let (step, progress, energy_prev) = update_step_adaptive(1.0, 2, 2.0, 1.0, &settings);
        let expected = 1.0 * settings.adaptive_decay;
        assert!((step - expected).abs() < 1e-9);
        assert_eq!(progress, 0);
        assert!((energy_prev - 2.0).abs() < 1e-12);
    }
}
