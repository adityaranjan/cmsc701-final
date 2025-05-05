use std::env;
use std::fs::File;
use std::io::Write;
use std::cmp::Ordering;
use minimizer_sa::shared::{MinimizerStringData, compute_minimizers, get_reference};


// TODO: add a terminal minimizer to the minimizer sequence?
// maybe not since a prefix of another sequence is less than that sequence with how rust does comparison


// buildsa function to work on minimizers (storing indices only)
fn buildsa(reference_path: &str, minimizer_k: usize, window_w: usize, output: &str) -> () {
    // Get the original reference sequence. This is needed for comparisons during sorting.
    let original_reference = get_reference(reference_path);

    if original_reference.is_empty() {
         eprintln!("Error: Reference sequence is empty.");
         return;
    }

    // Compute the minimizer sequence, which is now a vector of original genome positions.
    let minimizer_sequence = compute_minimizers(&original_reference, minimizer_k, window_w);

    if minimizer_sequence.is_empty() {
        eprintln!("Error: Could not compute minimizers. Check parameters or reference sequence length.");
        return;
    }

    // Build the suffix array on the *minimizer_sequence* (Vec<usize> of original positions).
    // The indices in minimizer_sa will point into the `minimizer_sequence` vector.
    let mut minimizer_sa: Vec<usize> = (0..minimizer_sequence.len()).collect();

    // Sort based on suffixes of the minimizer sequence.
    // The comparison now involves looking up k-mer strings in the original_reference
    // using the original genome positions stored in minimizer_sequence.
    minimizer_sa.sort_by(|idx1: &usize, idx2| {
        // Get the starting indices in the minimizer_sequence for the two suffixes being compared
        let mut sa_idx1_in_minimizers = *idx1;
        let mut sa_idx2_in_minimizers = *idx2;

        // Compare the suffixes element by element (where each element is a minimizer's k-mer)
        loop {
            // Check if we've reached the end of either suffix in the minimizer sequence
            let end_of_suffix1 = sa_idx1_in_minimizers >= minimizer_sequence.len();
            let end_of_suffix2 = sa_idx2_in_minimizers >= minimizer_sequence.len();

            match (end_of_suffix1, end_of_suffix2) {
                (true, true) => return Ordering::Equal, // Both suffixes end at the same time
                (true, false) => return Ordering::Less, // Suffix 1 is a prefix of Suffix 2
                (false, true) => return Ordering::Greater, // Suffix 2 is a prefix of Suffix 1
                (false, false) => {
                    // Get the original genome positions for the current minimizers
                    let pos1 = minimizer_sequence[sa_idx1_in_minimizers];
                    let pos2 = minimizer_sequence[sa_idx2_in_minimizers];

                    // Extract the k-mer strings from the original reference
                    // Ensure we don't go out of bounds when extracting k-mers
                    let kmer1 = &original_reference[pos1..(pos1 + minimizer_k).min(original_reference.len())];
                    let kmer2 = &original_reference[pos2..(pos2 + minimizer_k).min(original_reference.len())];


                    // Compare the k-mer strings
                    match kmer1.cmp(kmer2) {
                        Ordering::Equal => {
                            // K-mers are equal, continue to the next minimizer in the sequence
                            sa_idx1_in_minimizers += 1;
                            sa_idx2_in_minimizers += 1;
                        }
                        ordering => return ordering, // K-mers are different, return the comparison result
                    }
                }
            }
        }
    });

    // Store the new data structure
    let minimizer_string_data = MinimizerStringData {
        reference: original_reference, // Store original reference for querying
        minimizer_sequence, // Stores original genome positions
        minimizer_sa, // Indices into minimizer_sequence
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

    // Parse minimizer_k and window_w
    let minimizer_k: usize = args[2].parse().expect("minimizer_k must be an integer!");
    let window_w: usize = args[3].parse().expect("window_w must be an integer!");

    let output = &args[4];

    buildsa(reference_path, minimizer_k, window_w, output);
}
