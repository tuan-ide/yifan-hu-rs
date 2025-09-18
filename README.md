# Yifan Hu Graph Layout

This crate implements the multilevel spring-electrical graph layout described in
Yifan Hu's "Efficient and High Quality Force-Directed Graph Drawing". It
exposes both a reusable Rust library and a CLI binary for laying out undirected
graphs with straight-line edges.

## Using as a Library

Add the crate to your `Cargo.toml`:

```toml
[dependencies]
yifan-hu-graph-layout = { git = "https://github.com/tuan-ide/yifan-hu-graph-layout" }
```

Then call the multilevel layout from your code:

```rust
use yifan_hu_graph_layout::{Graph, LayoutSettings, multilevel_layout};

fn main() {
    // Build a simple line graph with four vertices.
    let mut graph = Graph::new(4);
    graph.add_edge(0, 1, 1.0);
    graph.add_edge(1, 2, 1.0);
    graph.add_edge(2, 3, 1.0);

    // Configure layout parameters. Most fields have sensible defaults.
    let mut settings = LayoutSettings::default();
    settings.max_iterations = 200;
    settings.tolerance = 1e-3;

    // Run the multilevel Yifan Hu layout.
    let result = multilevel_layout(&graph, &settings);

    for (idx, position) in result.positions.iter().enumerate() {
        println!("vertex {idx}: ({:.3}, {:.3})", position.x, position.y);
    }
}
```

## Using the CLI

```bash
cargo run -- --input graph.txt --output coords.csv
```

Graphs are expected as plain-text edge lists: the first non-comment line gives
the vertex count, followed by lines of `u v [weight]` with zero-based indices.
Coordinates are written as CSV.

## Development

- `cargo fmt`
- `cargo clippy --all-targets`
- `cargo test`
