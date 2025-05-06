use std::fs::File;
use std::cmp::Ordering;
use std::io::{BufRead, BufReader};
use serde::{Deserialize, Serialize};


#[derive(Serialize, Deserialize)]
pub struct MinimizerStringData {
    pub reference: String, // Store original reference for validation during querying
    pub minimizer_sequence: Vec<usize>, // stores minimizer indices in original genome sequence
    pub minimizer_sa: Vec<usize>, // suffix array for minimizer sequence
    pub minimizer_k: usize, // k for minimizers
    pub window_w: usize, // w for minimizers
}

// Function to compute the minimizer sequence (indices into original sequence)
// also collapses consecutive duplicate minimizers
pub fn compute_minimizers(sequence: &str, k: usize, w: usize) -> Vec<usize> {
    let mut minimizer_original_positions: Vec<usize> = Vec::new();

    let effective_len = sequence.len();

    // Store the original genome position of the last added minimizer
    let mut last_added_minimizer_original_pos: Option<usize> = None;

    // Slide a window of size w
    for i in 0..=(effective_len.saturating_sub(w)) {
        let window = &sequence[i..(std::cmp::min(i + w, effective_len))];

        let mut min_kmer_in_window: Option<&str> = None;
        let mut min_kmer_start_in_window = 0; // Start position of minimizer relative to window start

        // Find the lexicographically smallest k-mer in the current window
        for j in 0..=(window.len().saturating_sub(k)) {
            let current_kmer = &window[j..(j + k)];

            match min_kmer_in_window {
                None => {
                    min_kmer_in_window = Some(current_kmer);
                    min_kmer_start_in_window = j;
                }
                Some(existing_min) => {
                    if current_kmer < existing_min {
                        min_kmer_in_window = Some(current_kmer);
                        min_kmer_start_in_window = j;
                    }
                }
            }
        }

        if let Some(minimizer_kmer) = min_kmer_in_window {
            let current_minimizer_original_pos = i + min_kmer_start_in_window;

            let should_add = match last_added_minimizer_original_pos {
                None => true, // Always add the first minimizer found
                Some(last_pos) => {
                    // Extract the k-mer string for the last added minimizer
                    let last_kmer = &sequence[last_pos..(last_pos + k)];
                    // Collapse duplicates: add if the current minimizer k-mer is different from the last added k-mer
                    minimizer_kmer != last_kmer
                }
            };

            if should_add {
                minimizer_original_positions.push(current_minimizer_original_pos);
                last_added_minimizer_original_pos = Some(current_minimizer_original_pos);
            }
        }
    }

    minimizer_original_positions
}

// Helper function to compare two sequences of minimizer original positions
pub fn compare_minimizer_sequences(
    start_idx1: Option<usize>,
    start_idx2: Option<usize>,
    seq1_indices: &[usize],
    seq2_indices: &[usize],
    original_sequence1: &str,
    original_sequence2: &str,
    k: usize,
) -> Ordering {
    let mut idx1 = start_idx1.unwrap_or(0);
    let mut idx2 = start_idx2.unwrap_or(0);

    loop {
        // Check if we've reached the end of either sequence of indices
        let end_of_seq1 = idx1 >= seq1_indices.len();
        let end_of_seq2 = idx2 >= seq2_indices.len();

        match (end_of_seq1, end_of_seq2) {
            (true, true) => return Ordering::Equal, // Both sequences end at the same time
            (true, false) => return Ordering::Less, // Sequence 1 is a prefix of Sequence 2
            (false, true) => return Ordering::Greater, // Sequence 2 is a prefix of Sequence 1
            (false, false) => {
                // Get the original positions for the current minimizers
                let pos1 = seq1_indices[idx1];
                let pos2 = seq2_indices[idx2];

                // Extract the k-mer strings from their respective original sequences
                let kmer1 = &original_sequence1[pos1..(pos1 + k).min(original_sequence1.len())];
                let kmer2 = &original_sequence2[pos2..(pos2 + k).min(original_sequence2.len())];

                // Compare the k-mer strings
                match kmer1.cmp(kmer2) {
                    Ordering::Equal => {
                        idx1 += 1;
                        idx2 += 1;
                    }
                    ordering => return ordering,
                }
            }
        }
    }
}

pub fn get_reference(reference_path: &str) -> String {
    let file = File::open(reference_path).expect("file couldn't be opened!");
    let reader = BufReader::new(file);

    let result = reader
        .lines()
        .skip(1) // Skip the FASTA header line
        .map(|line| line.expect("failed to parse line!"))
        .collect::<String>();

    result
}
