use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use yifan_hu_rs::{Graph, LayoutSettings, Vec2, multilevel_layout};

fn main() -> Result<(), Box<dyn Error>> {
    let options = parse_args()?;
    let mut settings = LayoutSettings::default();
    if let Some(iterations) = options.max_iterations {
        settings.max_iterations = iterations;
    }
    if let Some(tolerance) = options.tolerance {
        settings.tolerance = tolerance;
    }
    if let Some(theta) = options.theta {
        settings.theta = theta;
    }
    if let Some(exponent) = options.repulsive_exponent {
        settings.repulsive_exponent = exponent;
    }
    if let Some(seed) = options.seed {
        settings.random_seed = Some(seed);
    }

    let graph = Graph::load_from_file(&options.input)?;
    let result = multilevel_layout(&graph, &settings);

    println!(
        "Laid out graph with {} vertices and {} edges in {} iterations.",
        graph.node_count(),
        graph.edge_count(),
        result.iterations
    );

    if let Some(output) = options.output {
        write_positions(&output, &result.positions)?;
        println!("Wrote coordinates to {}", output.display());
    } else {
        for (idx, pos) in result.positions.iter().enumerate() {
            println!("{idx} {:.6} {:.6}", pos.x, pos.y);
        }
    }

    Ok(())
}

struct CliOptions {
    input: PathBuf,
    output: Option<PathBuf>,
    max_iterations: Option<usize>,
    tolerance: Option<f64>,
    theta: Option<f64>,
    repulsive_exponent: Option<f64>,
    seed: Option<u64>,
}

fn parse_args() -> Result<CliOptions, Box<dyn Error>> {
    let mut args = env::args().skip(1);
    if env::args().any(|arg| arg == "--help" || arg == "-h") {
        print_usage();
        std::process::exit(0);
    }

    let mut input = None;
    let mut output = None;
    let mut max_iterations = None;
    let mut tolerance = None;
    let mut theta = None;
    let mut repulsive_exponent = None;
    let mut seed = None;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--input" => {
                let value = args.next().ok_or("Missing value for --input")?;
                input = Some(PathBuf::from(value));
            }
            "--output" => {
                let value = args.next().ok_or("Missing value for --output")?;
                output = Some(PathBuf::from(value));
            }
            "--iterations" => {
                let value = args.next().ok_or("Missing value for --iterations")?;
                max_iterations = Some(value.parse()?);
            }
            "--tolerance" => {
                let value = args.next().ok_or("Missing value for --tolerance")?;
                tolerance = Some(value.parse()?);
            }
            "--theta" => {
                let value = args.next().ok_or("Missing value for --theta")?;
                theta = Some(value.parse()?);
            }
            "--exponent" => {
                let value = args.next().ok_or("Missing value for --exponent")?;
                repulsive_exponent = Some(value.parse()?);
            }
            "--seed" => {
                let value = args.next().ok_or("Missing value for --seed")?;
                seed = Some(value.parse()?);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(&arg));
                } else if output.is_none() {
                    output = Some(PathBuf::from(arg));
                } else {
                    return Err(format!("Unrecognized argument: {arg}").into());
                }
            }
        }
    }

    let input = input.ok_or("No input graph specified")?;
    Ok(CliOptions {
        input,
        output,
        max_iterations,
        tolerance,
        theta,
        repulsive_exponent,
        seed,
    })
}

fn print_usage() {
    eprintln!(
        "Usage: yifan-hu-rs [--input path] [--output path] [--iterations n] [--tolerance tol] [--theta value] [--exponent value] [--seed value]"
    );
    eprintln!("You can also pass the input and output file paths positionally.");
    eprintln!(
        "Input format: first non-comment line contains number of vertices; subsequent lines contain edges as 'u v [weight]' with zero-based indices."
    );
}

fn write_positions(path: &PathBuf, positions: &[Vec2]) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "vertex,x,y")?;
    for (idx, pos) in positions.iter().enumerate() {
        writeln!(writer, "{idx},{:.6},{:.6}", pos.x, pos.y)?;
    }
    writer.flush()?;
    Ok(())
}
