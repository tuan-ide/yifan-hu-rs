use crate::graph::{Graph, Neighbor};
use rand::Rng;
use rand::seq::SliceRandom;
use std::collections::{HashMap, VecDeque};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum CoarseningStrategy {
    EdgeCollapse,
    Mivs,
    Hybrid,
}

#[derive(Debug, Clone)]
pub struct CoarsenResult {
    pub coarse: Graph,
    pub prolongation: Vec<Vec<(usize, f64)>>,
    pub used_strategy: CoarseningStrategy,
}

pub fn coarsen<R: Rng + ?Sized>(
    graph: &Graph,
    strategy: CoarseningStrategy,
    ratio_threshold: f64,
    rng: &mut R,
) -> Option<CoarsenResult> {
    if graph.node_count() <= 2 {
        return None;
    }
    match strategy {
        CoarseningStrategy::EdgeCollapse => edge_collapse(graph, rng),
        CoarseningStrategy::Mivs => mivs(graph, rng),
        CoarseningStrategy::Hybrid => {
            let ec_result = edge_collapse(graph, rng);
            if let Some(result) = ec_result {
                let ratio = result.coarse.node_count() as f64 / graph.node_count() as f64;
                if ratio <= ratio_threshold {
                    return Some(result);
                }
            }
            mivs(graph, rng)
        }
    }
}

fn edge_collapse<R: Rng + ?Sized>(graph: &Graph, rng: &mut R) -> Option<CoarsenResult> {
    let n = graph.node_count();
    if n <= 1 {
        return None;
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.shuffle(rng);
    let mut matched = vec![false; n];
    let mut coarse_nodes: Vec<Vec<usize>> = Vec::new();
    let mut fine_to_coarse = vec![usize::MAX; n];

    for &u in &order {
        if matched[u] {
            continue;
        }
        matched[u] = true;
        let mut best_neighbor: Option<usize> = None;
        let mut best_weight = -f64::INFINITY;
        let mut neighbor_order: Vec<&Neighbor> = graph.neighbors(u).iter().collect();
        neighbor_order.shuffle(rng);
        for neigh in neighbor_order {
            if matched[neigh.target] {
                continue;
            }
            // Prefer heavier edges; tie broken by randomness via shuffle.
            if neigh.weight > best_weight {
                best_weight = neigh.weight;
                best_neighbor = Some(neigh.target);
            }
        }

        let coarse_idx = coarse_nodes.len();
        match best_neighbor {
            Some(v) => {
                matched[v] = true;
                coarse_nodes.push(vec![u, v]);
                fine_to_coarse[u] = coarse_idx;
                fine_to_coarse[v] = coarse_idx;
            }
            None => {
                coarse_nodes.push(vec![u]);
                fine_to_coarse[u] = coarse_idx;
            }
        }
    }

    let coarse_count = coarse_nodes.len();
    if coarse_count == n {
        // No effective coarsening happened.
        return None;
    }

    let mut prolongation = vec![Vec::<(usize, f64)>::new(); n];
    for (coarse_idx, fine_group) in coarse_nodes.iter().enumerate() {
        for &fine_idx in fine_group {
            prolongation[fine_idx].push((coarse_idx, 1.0));
        }
    }

    let coarse = build_coarse_graph(graph, coarse_count, &prolongation);
    Some(CoarsenResult {
        coarse,
        prolongation,
        used_strategy: CoarseningStrategy::EdgeCollapse,
    })
}

fn mivs<R: Rng + ?Sized>(graph: &Graph, rng: &mut R) -> Option<CoarsenResult> {
    let n = graph.node_count();
    if n <= 1 {
        return None;
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.shuffle(rng);

    let mut in_set = vec![false; n];
    let mut blocked = vec![false; n];

    for &v in &order {
        if blocked[v] {
            continue;
        }
        in_set[v] = true;
        blocked[v] = true;
        for neighbor in graph.neighbors(v) {
            blocked[neighbor.target] = true;
        }
    }

    let set_indices: Vec<usize> = (0..n).filter(|&i| in_set[i]).collect();
    if set_indices.is_empty() {
        return None;
    }
    let coarse_count = set_indices.len();
    if coarse_count == n {
        return None;
    }

    let mut index_map = vec![None; n];
    for (idx, &v) in set_indices.iter().enumerate() {
        index_map[v] = Some(idx);
    }

    let mut prolongation = vec![Vec::<(usize, f64)>::new(); n];

    // Assign members of the independent set directly.
    for (coarse_idx, &fine_idx) in set_indices.iter().enumerate() {
        prolongation[fine_idx].push((coarse_idx, 1.0));
    }

    // Assign remaining vertices by averaging neighbors in the set.
    for v in 0..n {
        if in_set[v] {
            continue;
        }
        let mut candidates: Vec<usize> = graph
            .neighbors(v)
            .iter()
            .filter_map(|neigh| index_map[neigh.target])
            .collect();

        if candidates.is_empty() {
            // Fall back to BFS within distance 3 to find nearest independent nodes.
            candidates = nearest_mivs_nodes(graph, v, &index_map, 3);
        }

        if candidates.is_empty() {
            // As a last resort, choose a random representative.
            let choose = order[0];
            if let Some(idx) = index_map[choose] {
                candidates.push(idx);
            } else if let Some(idx) = index_map.iter().flatten().next().copied() {
                candidates.push(idx);
            } else {
                return None;
            }
        }

        candidates.sort_unstable();
        candidates.dedup();
        let weight = 1.0 / candidates.len() as f64;
        for coarse_idx in candidates {
            prolongation[v].push((coarse_idx, weight));
        }
    }

    let coarse = build_coarse_graph(graph, coarse_count, &prolongation);
    Some(CoarsenResult {
        coarse,
        prolongation,
        used_strategy: CoarseningStrategy::Mivs,
    })
}

fn nearest_mivs_nodes(
    graph: &Graph,
    start: usize,
    index_map: &Vec<Option<usize>>,
    max_depth: usize,
) -> Vec<usize> {
    let mut visited = vec![false; graph.node_count()];
    let mut queue = VecDeque::new();
    queue.push_back((start, 0usize));
    visited[start] = true;
    let mut found = Vec::new();
    let mut best_depth = None;

    while let Some((node, depth)) = queue.pop_front() {
        if depth > max_depth {
            break;
        }
        if let Some(idx) = index_map[node] {
            if best_depth.map_or(true, |best| depth <= best) {
                if best_depth.is_none() {
                    best_depth = Some(depth);
                }
                found.push(idx);
                continue;
            }
        }
        if Some(depth) == best_depth {
            continue;
        }
        for neighbor in graph.neighbors(node) {
            if !visited[neighbor.target] {
                visited[neighbor.target] = true;
                queue.push_back((neighbor.target, depth + 1));
            }
        }
    }
    found
}

fn build_coarse_graph(
    graph: &Graph,
    coarse_count: usize,
    prolongation: &[Vec<(usize, f64)>],
) -> Graph {
    let mut coarse = Graph::new(coarse_count);
    let mut coarse_weights = vec![0.0; coarse_count];
    let mut edge_weights: HashMap<(usize, usize), f64> = HashMap::new();

    for (fine_idx, weights) in prolongation.iter().enumerate() {
        let fine_weight = graph.node_weight(fine_idx);
        for &(coarse_idx, w) in weights {
            coarse_weights[coarse_idx] += fine_weight * w;
        }
    }

    for (idx, weight) in coarse_weights.iter().enumerate() {
        coarse.set_node_weight(idx, *weight);
    }

    for u in 0..graph.node_count() {
        for neighbor in graph.neighbors(u) {
            if neighbor.target < u {
                continue;
            }
            let fine_v = neighbor.target;
            let weight = neighbor.weight;
            for &(cu, wu) in &prolongation[u] {
                for &(cv, wv) in &prolongation[fine_v] {
                    if cu == cv {
                        continue;
                    }
                    let (a, b) = if cu < cv { (cu, cv) } else { (cv, cu) };
                    let contribution = weight * wu * wv;
                    *edge_weights.entry((a, b)).or_insert(0.0) += contribution;
                }
            }
        }
    }

    for ((u, v), weight) in edge_weights {
        coarse.add_edge(u, v, weight.max(1e-9));
    }

    coarse
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn path_graph(n: usize) -> Graph {
        let mut g = Graph::new(n);
        for i in 0..n - 1 {
            g.add_edge(i, i + 1, 1.0);
        }
        g
    }

    #[test]
    fn edge_collapse_reduces_vertex_count() {
        let graph = path_graph(6);
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let result =
            coarsen(&graph, CoarseningStrategy::EdgeCollapse, 0.75, &mut rng).expect("coarsened");
        assert!(result.coarse.node_count() < graph.node_count());
        assert_eq!(result.used_strategy, CoarseningStrategy::EdgeCollapse);
    }

    #[test]
    fn hybrid_falls_back_to_mivs_when_ratio_exceeds_threshold() {
        let graph = path_graph(8);
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let result = coarsen(&graph, CoarseningStrategy::Hybrid, 0.1, &mut rng).expect("coarsened");
        assert_eq!(result.used_strategy, CoarseningStrategy::Mivs);
        assert!(result.coarse.node_count() < graph.node_count());
    }
}
