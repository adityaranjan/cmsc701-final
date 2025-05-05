use std::env;
use std::fs::File;
use std::io::Write;
use std::cmp::Ordering;
use minimizer_sa::shared::{MinimizerStringData, get_reference};


// TODO: add a terminal minimizer to the minimizer sequence?
// maybe not since a prefix of another sequence is less than that sequence with rust comparison


// Function to compute the minimizer sequence and map to original positions
// Computes lexicographical minimizers
fn compute_minimizers(sequence: &str, k: usize, w: usize) -> (Vec<String>, Vec<usize>) {
    let mut minimizer_seq: Vec<String> = Vec::new();
    let mut original_pos: Vec<usize> = Vec::new();

    // Process the sequence
    let effective_len = sequence.len();

    if effective_len < k {
        // Sequence is too short to even form a k-mer
        return (minimizer_seq, original_pos);
    }

    // Slide a window of size w
    for i in 0..=(effective_len.saturating_sub(w)) {
        let window = &sequence[i..std::cmp::min(i + w, effective_len)];

        if window.len() < k {
            // Window is too short to contain a k-mer of size k
            break; // Should not happen if loop bounds are correct
        }

        let mut min_kmer: Option<&str> = None;
        let mut min_kmer_start_in_window = 0; // Start position of minimizer relative to window start

        // Find the lexicographically smallest k-mer in the current window
        for j in 0..=(window.len().saturating_sub(k)) {
            let current_kmer = &window[j..j + k];

            match min_kmer {
                None => {
                    // First k-mer in the window becomes the initial minimum
                    min_kmer = Some(current_kmer);
                    min_kmer_start_in_window = j;
                }
                Some(existing_min) => {
                    // Compare current k-mer with the current minimum
                    if current_kmer.cmp(existing_min) < Ordering::Less {
                        min_kmer = Some(current_kmer);
                        min_kmer_start_in_window = j;
                    }
                }
            }
        }

        // If a valid minimizer was found in the window (it should be if w >= k)
        if let Some(minimizer) = min_kmer {
            minimizer_seq.push(minimizer.to_string());
            // The original genome position for this minimizer is the start of the window (i)
            // plus the start of the minimizer within that window.
            original_pos.push(i + min_kmer_start_in_window);
        }
    }

    (minimizer_seq, original_pos)
}


// buildsa function to work on minimizers
fn buildsa(reference_path: &str, minimizer_k: usize, window_w: usize, output: &str) -> () {
    let original_reference = get_reference(reference_path);

    // Compute the minimizer sequence and the mapping
    let (minimizer_sequence, minimizer_to_genome_pos) =
        compute_minimizers(&original_reference, minimizer_k, window_w);

    // Check if minimizer sequence was generated
    if minimizer_sequence.is_empty() {
        eprintln!("Error: Could not compute minimizers. Check parameters or reference sequence.");
        // Consider if you want to save an empty file or indicate failure differently
        return;
    }

    // Build the suffix array on the *minimizer sequence*
    // The indices in minimizer_sa will point into the `minimizer_sequence` vector.
    let mut minimizer_sa: Vec<usize> = (0..minimizer_sequence.len()).collect();

    // Sort based on suffixes of the minimizer sequence
    minimizer_sa.sort_by(|idx1, idx2| {
        // Compare the suffixes of the minimizer_sequence starting at *idx1* and *idx2*
        (&minimizer_sequence[*idx1..]).cmp(&minimizer_sequence[*idx2..])
    });

    // Store the new data structure
    let minimizer_string_data = MinimizerStringData {
        reference: original_reference,
        minimizer_sequence,
        minimizer_to_genome_pos,
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
    println!("Minimizer sequence length: {}", minimizer_string_data.minimizer_sequence.len());
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
