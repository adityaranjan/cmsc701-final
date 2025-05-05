use std::env;
use std::fs::File;
use std::cmp::Ordering;
use std::io::{BufRead, BufReader, Read, Write};
use minimizer_sa::shared::{MinimizerStringData, compute_minimizers};


fn get_data(index: &str) -> MinimizerStringData {
    let mut file = File::open(index).expect("failed to open the index file!");
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)
        .expect("failed to read the index file!");

    bincode::deserialize(&buffer).expect("failed to deserialize the index data!")
}


fn create_terminal_minimizer_pattern(k: usize) -> Vec<String> {
    // A single k-mer string of '{' characters, which is lexicographically
    // greater than any k-mer composed of A, C, G, T.
    let terminal_kmer: String = "{".repeat(k);
    vec![terminal_kmer]
}


// Helper function to compare two sequences of minimizer original positions
// by extracting and comparing the actual k-mer strings from their respective original sequences.
fn compare_minimizer_sequences(
    seq1_indices: &[usize], // Slice of original positions for sequence 1 (e.g., reference suffix)
    seq2_indices: &[usize], // Slice of original positions for sequence 2 (e.g., query pattern)
    original_sequence1: &str, // The original DNA sequence for sequence 1 (e.g., original reference)
    original_sequence2: &str, // The original DNA sequence for sequence 2 (e.g., original query)
    k: usize, // Minimizer k-mer size
) -> Ordering {
    let mut idx1 = 0;
    let mut idx2 = 0;

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
                // Ensure we don't go out of bounds when extracting k-mers
                let kmer1 = &original_sequence1[pos1..(pos1 + k).min(original_sequence1.len())];
                let kmer2 = &original_sequence2[pos2..(pos2 + k).min(original_sequence2.len())];

                // Compare the k-mer strings
                match kmer1.cmp(kmer2) {
                    Ordering::Equal => {
                        // K-mers are equal, continue to the next minimizer index
                        idx1 += 1;
                        idx2 += 1;
                    }
                    ordering => return ordering, // K-mers are different, return the comparison result
                }
            }
        }
    }
}


// Modified bin_search to work on sequences of original positions.
// This function finds the lower bound index in the minimizer_sa for a given pattern (sequence of original positions).
fn bin_search(
    left: usize,
    right: usize,
    minimizer_sequence_indices: &Vec<usize>, // The reference minimizer sequence (original positions)
    minimizer_sa: &Vec<usize>,      // The suffix array on minimizer_sequence_indices
    pattern_minimizer_indices: &Vec<usize>, // The query minimizer sequence (original positions in query)
    original_reference: &str, // Original reference sequence for k-mer extraction
    original_query_sequence: &str, // Original query sequence for k-mer extraction
    k: usize, // Minimizer k-mer size
) -> usize {
    let mut l = left;
    let mut r = right;

    while l < r {
        let m = l + (r - l) / 2; // Use this to prevent potential overflow

        // Get the slice of original positions for the suffix in the reference minimizer sequence
        let suffix_in_minimizers_indices = &minimizer_sequence_indices[minimizer_sa[m]..];

        // Compare the suffix (sequence of indices from reference) with the pattern (sequence of indices from query)
        // using the helper function to extract and compare k-mers from their respective original sequences.
        match compare_minimizer_sequences(
            suffix_in_minimizers_indices, // Indices from the reference minimizer sequence
            pattern_minimizer_indices, // Indices from the query minimizer sequence
            original_reference, // Use original reference for suffix k-mers
            original_query_sequence, // Use original query for pattern k-mers
            k,
        ) {
            Ordering::Less => {
                // Suffix is less than pattern, need to search in the right half
                l = m + 1;
            }
            Ordering::Equal => {
                // Found a match, but it might not be the first. Search in the left half for the lower bound.
                r = m;
            }
            Ordering::Greater => {
                // Suffix is greater than pattern, need to search in the left half
                r = m;
            }
        }
    }

    // l is now the lower bound (first index in minimizer_sa where the suffix is >= pattern)
    l
}


// Modified process_query to handle minimizer-space search and report all potential matches
// Removed validation and false positive rate calculation.
fn process_query(
    data: &MinimizerStringData,
    original_query_sequence: String, // The original query DNA sequence
    query_name: &str,
    output_file: &mut File,
) -> () {
    // Transform the original query sequence into its minimizer sequence (indices into query)
    let query_minimizer_indices = compute_minimizers(
        &original_query_sequence,
        data.minimizer_k,
        data.window_w,
    );

    // If the query minimizer sequence is empty, it cannot match anything.
    if query_minimizer_indices.is_empty() {
         let output_string = format!("{}\t0", query_name); // Report 0 occurrences
         writeln!(output_file, "{}", output_string).expect("failed to write to the output file!");
         return;
    }

    // --- Find the range [l, r) in the minimizer_sa using two binary searches ---

    // 1. Find the lower bound (l): First suffix >= query_minimizer_indices
    let l = bin_search(
        0, // Search the entire minimizer_sa
        data.minimizer_sa.len(),
        &data.minimizer_sequence, // Reference minimizer indices
        &data.minimizer_sa,
        &query_minimizer_indices, // Query minimizer indices
        &data.reference, // Original reference for comparisons
        &original_query_sequence, // Original query for comparisons
        data.minimizer_k,
    );

    // 2. Find the upper bound (r): First suffix > query_minimizer_indices
    // We need to find the first suffix in the minimizer SA that is lexicographically
    // strictly greater than the query minimizer sequence.
    // This is equivalent to finding the lower bound of the "next" sequence after the query.
    // Since we are comparing sequences of indices by looking up k-mers, the simplest
    // way to find the upper bound is to iterate from 'l' and stop when the suffix
    // no longer starts with the query minimizer sequence.

    let mut r = l;
    while r < data.minimizer_sa.len() {
        // Get the slice of original positions for the suffix in the reference minimizer sequence
        let suffix_in_minimizers_indices = &data.minimizer_sequence[data.minimizer_sa[r]..];

        // Check if the suffix starts with the query minimizer sequence using k-mer comparison
        // We compare the suffix slice with the full query minimizer index sequence.
        // If the comparison is Ordering::Equal or Ordering::Greater, it means the suffix
        // starts with or is equal to the query minimizer sequence. We want the first
        // index where it is strictly greater.
        // A simpler check is to see if the prefix of the suffix matches the query.
        if suffix_in_minimizers_indices.len() >= query_minimizer_indices.len() &&
           compare_minimizer_sequences(
               &suffix_in_minimizers_indices[0..query_minimizer_indices.len()], // Compare prefix of suffix
               &query_minimizer_indices, // with the full query sequence
               &data.reference, // Use original reference for suffix k-mers
               &original_query_sequence, // Use original query for pattern k-mers
               data.minimizer_k
           ) == Ordering::Equal
        {
            r += 1; // The suffix starts with the query, continue
        } else {
            break; // Suffix no longer starts with the query
        }
    }


    // --- Report potential matches found in minimizer space ---

    let potential_match_positions: Vec<usize> = (l..r)
        .map(|i| data.minimizer_sequence[data.minimizer_sa[i]]) // Get the original genome position
        .collect();

    // Output the results
    let mut output_string = String::new();
    output_string.push_str(query_name);
    output_string.push_str("\t");
    output_string.push_str(&potential_match_positions.len().to_string()); // Report the count of potential matches

    for pos in potential_match_positions {
        output_string.push_str("\t");
        output_string.push_str(&pos.to_string()); // Report the original genome positions
    }

    writeln!(output_file, "{}", output_string).expect("failed to write to the output file!");
}

// Modified querysa function
fn querysa(index: &str, queries: &str, output: &str) -> () {
    // get the serialized reference data
    let data = get_data(index);

    // open the queries file
    let file = File::open(queries).expect("queries file couldn't be opened!!");
    let reader = BufReader::new(file);

    let mut output_file = File::create(output).expect("failed to create the file!");

    let mut num_records = 0;
    let mut curr_sequence: String;
    let mut curr_query = String::new();
    let mut curr_sequence_vec: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line.expect("failed to read the line!");

        if let Some(first_char) = line.chars().next() {
            if first_char == '>' {
                if num_records > 0 {
                    curr_sequence = curr_sequence_vec.join("");
                    curr_sequence_vec.clear();
                    // Call process_query
                    process_query(
                        &data,
                        curr_sequence,
                        &curr_query,
                        &mut output_file,
                    );
                }

                curr_query = (&line[1..]).to_string();
                num_records += 1;
            } else {
                curr_sequence_vec.push(line);
            }
        }
    }

    // Process the last query
    curr_sequence = curr_sequence_vec.join("");
    curr_sequence_vec.clear();
    process_query(
        &data,
        curr_sequence,
        &curr_query,
        &mut output_file,
    );
}
fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 4 {
        eprintln!("Usage: {} <index_path> <queries_path> <output_path>", args[0]);
        return;
    }

    let index = &args[1];
    let queries = &args[2];
    let output = &args[3];

    querysa(index, queries, output);
}
