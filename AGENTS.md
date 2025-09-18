# Repository Guidelines

## Project Structure & Module Organization
- `Cargo.toml` defines the Rust 2024 crate and external dependencies; keep metadata and edition settings up to date.
- `src/main.rs` hosts the binary entry point; introduce helper modules under `src/` when the layout algorithm grows, and expose shared logic via `mod` declarations.
- `graph_draw_small.pdf` is a reference asset for expected graph aesthetics—use it to validate visual output when tuning forces or cooling schedules.

## Build, Test, and Run Commands
- `cargo build` compiles the crate in debug mode; run before pushing to catch compiler errors locally.
- `cargo run --release -- <args>` executes the layout binary with optimized settings; reserve for timing-sensitive experiments.
- `cargo test` executes unit and integration tests.
- `cargo fmt` and `cargo clippy --all-targets` keep formatting consistent and surface lints before review.

## Coding Style & Naming Conventions
- Stick with `rustfmt` defaults (4-space indentation, trailing commas where allowed); never hand-format.
- Use `snake_case` for functions and variables, `UpperCamelCase` for types, and uppercase `CONSTANTS` per Rust idioms.
- Favor iterator-driven code and expressive structs over deeply nested loops to model forces and constraints.

## Testing Guidelines
- Place fast unit tests inside `#[cfg(test)]` modules near the code they cover; integration tests belong under `tests/`.
- When adding algorithms, craft deterministic fixtures (e.g., small adjacency lists) and assert on computed coordinates or energy deltas.
- Aim for meaningful coverage of convergence steps and panic paths before requesting review.

## Commit & Pull Request Guidelines
- Mirror the existing `Init project` style: concise, imperative subjects under 60 characters (e.g., `Add cooling schedule tuner`).
- Reference related issues in the body (`Refs #12`) and call out user-visible changes or new CLI flags.
- PRs should summarize algorithmic changes, include profiling notes when performance shifts, and attach updated visuals if `graph_draw_small.pdf` expectations change.
