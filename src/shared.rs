use std::fs::File;
use std::cmp::Ordering;
use std::io::{BufRead, BufReader};
use serde::{Deserialize, Serialize};


// Define the struct for minimizer-space SA data
#[derive(Serialize, Deserialize)]
pub struct MinimizerStringData {
    pub reference: String, // Store original reference for validation during querying
    // minimizer_sequence now stores indices into the original genome sequence
    // where the selected minimizers start.
    pub minimizer_sequence: Vec<usize>,
    // minimizer_to_genome_pos is no longer needed as minimizer_sequence *is* the mapping
    // pub minimizer_to_genome_pos: Vec<usize>,
    // minimizer_sa stores indices into the *minimizer_sequence* (Vec<usize>)
    pub minimizer_sa: Vec<usize>,
    pub minimizer_k: usize, // k for minimizers
    pub window_w: usize, // w for minimizers
}

// Function to compute the minimizer sequence (indices into original sequence)
// and collapse consecutive duplicates based on their k-mer sequence.
// Returns a vector of original genome start positions of the selected minimizers.
pub fn compute_minimizers(sequence: &str, k: usize, w: usize) -> Vec<usize> {
    let mut minimizer_original_positions: Vec<usize> = Vec::new();

    let effective_len = sequence.len();

    if effective_len < k {
        return minimizer_original_positions;
    }

    // Store the original genome position of the last added minimizer
    let mut last_added_minimizer_original_pos: Option<usize> = None;

    // Slide a window of size w
    for i in 0..=(effective_len.saturating_sub(w)) {
        let window = &sequence[i..(std::cmp::min(i + w, effective_len))];

        if window.len() < k {
            break;
        }

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
                    if current_kmer.cmp(existing_min) < Ordering::Less {
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
                    // Add if the current minimizer k-mer is different from the last added k-mer
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

pub fn lcp(s1: &str, s2: &str) -> f64 {
    let mut i: usize = 0;

    while i < s1.len() && i < s2.len() && &s1[i..i + 1] == &s2[i..i + 1] {
        i += 1;
    }

    i as f64
}
