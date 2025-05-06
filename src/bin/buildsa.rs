use std::env;
use std::fs::File;
use std::io::Write;
use minimizer_sa::shared::{MinimizerStringData, compute_minimizers, compare_minimizer_sequences, get_reference};


// buildsa function to work on minimizers
fn buildsa(reference_path: &str, minimizer_k: usize, window_w: usize, output: &str) -> () {
    let mut original_reference = get_reference(reference_path);

    // Compute the minimizer sequence, which is now a vector of original genome positions.
    let mut minimizer_sequence = compute_minimizers(&original_reference, minimizer_k, window_w);

    // add terminal k-mer to original_reference and correspnding index to minimizer_sequence
    original_reference.push_str(&"$".repeat(minimizer_k));
    minimizer_sequence.push(original_reference.len() - minimizer_k);

    // Build the suffix array on the minimizer_sequence (Vec<usize> of original positions).
    // The indices in minimizer_sa will point into the minimizer_sequence vector.
    let mut minimizer_sa: Vec<usize> = (0..minimizer_sequence.len()).collect();

    // Sort based on suffixes of the minimizer sequence.
    // The comparison involves looking up k-mer strings in the original_reference
    // using the original genome positions stored in minimizer_sequence.
    minimizer_sa.sort_by(|idx1: &usize, idx2| {
        compare_minimizer_sequences(
            Some(*idx1),
            Some(*idx2),
             &minimizer_sequence[..], 
             &minimizer_sequence[..],
             &original_reference, 
             &original_reference, 
             minimizer_k,
        )
    });

    let minimizer_string_data = MinimizerStringData {
        reference: original_reference,
        minimizer_sequence,
        minimizer_sa,
        minimizer_k,
        window_w,
    };

    // Serialize and save
    let mut output_file = File::create(output).expect("output file creation failed!");
    let serialized_data = bincode::serialize(&minimizer_string_data).expect("failed to serialize!");
    output_file
        .write_all(&serialized_data)
        .expect("writing the data to the output file failed!");

    println!("Minimizer-space suffix array built successfully!");
    println!("Original sequence length: {}", minimizer_string_data.reference.len());
    println!("Minimizer sequence length (indices stored): {}", minimizer_string_data.minimizer_sequence.len());
    println!("Minimizer k: {}", minimizer_k);
    println!("Window w: {}", window_w);
}
fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: {} <reference_path> <minimizer_k> <window_w> <output_path>", args[0]);
        return;
    }

    let reference_path = &args[1];

    let minimizer_k: usize = args[2].parse().expect("minimizer_k must be an integer!");
    let window_w: usize = args[3].parse().expect("window_w must be an integer!");

    let output = &args[4];

    buildsa(reference_path, minimizer_k, window_w, output);
}
