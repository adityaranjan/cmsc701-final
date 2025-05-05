use std::fs::File;
use std::io::{BufRead, BufReader};
use serde::{Deserialize, Serialize};


// Define the struct for minimizer-space SA data
#[derive(Serialize, Deserialize)]
pub struct MinimizerStringData {
    pub reference: String, // Original genome sequence
    pub minimizer_sequence: Vec<String>, // Sequence of minimizer k-mers
    pub minimizer_to_genome_pos: Vec<usize>, // Mapping: index in minimizer_sequence -> start pos in original genome
    pub minimizer_sa: Vec<usize>, // Suffix array on minimizer_sequence (indices are into minimizer_sequence)
    pub minimizer_k: usize, // k for minimizers
    pub window_w: usize, // w for minimizers
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
