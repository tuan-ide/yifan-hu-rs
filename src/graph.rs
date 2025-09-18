use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug)]
pub enum GraphError {
    Io(std::io::Error),
    Parse(String),
}

impl From<std::io::Error> for GraphError {
    fn from(err: std::io::Error) -> Self {
        GraphError::Io(err)
    }
}

impl std::fmt::Display for GraphError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GraphError::Io(err) => write!(f, "IO error: {}", err),
            GraphError::Parse(msg) => write!(f, "Parse error: {}", msg),
        }
    }
}

impl std::error::Error for GraphError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            GraphError::Io(err) => Some(err),
            GraphError::Parse(_) => None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Neighbor {
    pub target: usize,
    pub weight: f64,
}

#[derive(Clone, Debug)]
pub struct Node {
    pub weight: f64,
    pub neighbors: Vec<Neighbor>,
}

#[derive(Clone, Debug)]
pub struct Graph {
    pub nodes: Vec<Node>,
}

impl Graph {
    pub fn new(node_count: usize) -> Self {
        Self {
            nodes: (0..node_count)
                .map(|_| Node {
                    weight: 1.0,
                    neighbors: Vec::new(),
                })
                .collect(),
        }
    }

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    pub fn edge_count(&self) -> usize {
        self.nodes.iter().map(|n| n.neighbors.len()).sum::<usize>() / 2
    }

    pub fn add_edge(&mut self, u: usize, v: usize, weight: f64) {
        assert!(
            u < self.nodes.len() && v < self.nodes.len(),
            "vertex index out of bounds"
        );
        if u == v {
            return;
        }
        if !self.nodes[u].neighbors.iter().any(|n| n.target == v) {
            self.nodes[u].neighbors.push(Neighbor { target: v, weight });
        }
        if !self.nodes[v].neighbors.iter().any(|n| n.target == u) {
            self.nodes[v].neighbors.push(Neighbor { target: u, weight });
        }
    }

    pub fn load_from_file(path: &Path) -> Result<Self, GraphError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let mut node_count = None;
        let mut graph: Option<Graph> = None;

        for (line_idx, line) in reader.lines().enumerate() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if node_count.is_none() {
                if parts.len() < 1 {
                    return Err(GraphError::Parse(format!(
                        "Expected node count on line {}",
                        line_idx + 1
                    )));
                }
                let n: usize = parts[0].parse().map_err(|_| {
                    GraphError::Parse(format!(
                        "Invalid node count '{}' on line {}",
                        parts[0],
                        line_idx + 1
                    ))
                })?;
                if n == 0 {
                    return Err(GraphError::Parse(
                        "Graph must have at least one vertex".into(),
                    ));
                }
                graph = Some(Graph::new(n));
                node_count = Some(n);
                continue;
            }

            let g = graph.as_mut().expect("graph initialized");
            if parts.len() < 2 {
                return Err(GraphError::Parse(format!(
                    "Expected edge definition on line {}",
                    line_idx + 1
                )));
            }

            let u: usize = parts[0].parse().map_err(|_| {
                GraphError::Parse(format!(
                    "Invalid vertex '{}' on line {}",
                    parts[0],
                    line_idx + 1
                ))
            })?;
            let v: usize = parts[1].parse().map_err(|_| {
                GraphError::Parse(format!(
                    "Invalid vertex '{}' on line {}",
                    parts[1],
                    line_idx + 1
                ))
            })?;

            let weight = if parts.len() >= 3 {
                parts[2].parse().map_err(|_| {
                    GraphError::Parse(format!(
                        "Invalid weight '{}' on line {}",
                        parts[2],
                        line_idx + 1
                    ))
                })?
            } else {
                1.0
            };

            if u >= g.node_count() || v >= g.node_count() {
                return Err(GraphError::Parse(format!(
                    "Vertex index out of bounds on line {}",
                    line_idx + 1
                )));
            }
            g.add_edge(u, v, weight);
        }

        graph.ok_or_else(|| GraphError::Parse("File did not contain a node count".into()))
    }

    pub fn neighbors(&self, node: usize) -> &Vec<Neighbor> {
        &self.nodes[node].neighbors
    }

    pub fn node_weight(&self, node: usize) -> f64 {
        self.nodes[node].weight
    }

    pub fn set_node_weight(&mut self, node: usize, weight: f64) {
        self.nodes[node].weight = weight;
    }

    pub fn approximate_diameter(&self) -> usize {
        if self.node_count() <= 1 {
            return 0;
        }
        let mut start = 0;
        let mut best_dist = 0;

        loop {
            let (target, dist) = self.bfs_farthest(start);
            if dist <= best_dist {
                return best_dist;
            }
            start = target;
            best_dist = dist;
        }
    }

    fn bfs_farthest(&self, start: usize) -> (usize, usize) {
        let mut visited = vec![false; self.node_count()];
        let mut queue = VecDeque::new();
        queue.push_back((start, 0usize));
        visited[start] = true;
        let mut farthest = (start, 0usize);

        while let Some((node, dist)) = queue.pop_front() {
            if dist > farthest.1 {
                farthest = (node, dist);
            }
            for neighbor in &self.nodes[node].neighbors {
                if !visited[neighbor.target] {
                    visited[neighbor.target] = true;
                    queue.push_back((neighbor.target, dist + 1));
                }
            }
        }

        farthest
    }
}
