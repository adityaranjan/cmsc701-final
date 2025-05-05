use std::env;
use std::fs::File;
use std::cmp::Ordering;
use std::io::{BufRead, BufReader, Read, Write};
use minimizer_sa::shared::MinimizerStringData;


// get_data to load MinimizerStringData
fn get_data(index: &str) -> MinimizerStringData {
    let mut file = File::open(index).expect("failed to open the index file!");
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)
        .expect("failed to read the index file!");

    bincode::deserialize(&buffer).expect("failed to deserialize the index data!")
}

// Function to compute the minimizer sequence for a query pattern
// (Similar to the builder's compute_minimizers, but only returns the sequence)
fn compute_minimizers_for_query(sequence: &str, k: usize, w: usize) -> Vec<String> {
    let mut minimizer_seq = Vec::new();

    // Process the sequence
    let effective_len = sequence.len();

    if effective_len < k {
        return minimizer_seq; // Sequence too short to form a k-mer
    }

    for i in 0..=(effective_len.saturating_sub(w)) {
        let window = &sequence[i..std::cmp::min(i + w, effective_len)];

         if window.len() < k {
            break; // Window too short for k-mer
        }

        let mut min_kmer: Option<&str> = None;

        for j in 0..=(window.len().saturating_sub(k)) {
            let current_kmer = &window[j..j + k];

            match min_kmer {
                None => {
                    min_kmer = Some(current_kmer);
                }
                Some(existing_min) => {
                    if current_kmer.cmp(existing_min) < Ordering::Less {
                        min_kmer = Some(current_kmer);
                    }
                }
            }
        }

        if let Some(minimizer) = min_kmer {
            minimizer_seq.push(minimizer.to_string());
        }
    }

    minimizer_seq
}

// Helper function to create a "terminal minimizer" sequence for upper bound search
fn create_terminal_minimizer_pattern(k: usize) -> Vec<String> {
    // A single k-mer string of '{' characters, which is lexicographically
    // greater than any k-mer composed of A, C, G, T.
    let terminal_kmer = "{".repeat(k);
    vec![terminal_kmer]
}


// bin_search to work on the minimizer sequence
// This function finds the lower bound index in the minimizer_sa for a given pattern (sequence of minimizers)
fn bin_search(
    left: usize,
    right: usize,
    minimizer_sequence: &Vec<String>, // The reference minimizer sequence
    minimizer_sa: &Vec<usize>,      // The suffix array on minimizers
    pattern_minimizers: &Vec<String>, // The query minimizer sequence (pattern)
) -> usize {
    let mut l = left;
    let mut r = right;

    while l < r {
        let m = l + (r - l) / 2; // Use this to prevent potential overflow

        // Compare the suffix of the minimizer_sequence starting at minimizer_sa[m]
        // with the pattern_minimizers sequence.
        let suffix_in_minimizers = &minimizer_sequence[minimizer_sa[m]..];

        match suffix_in_minimizers.cmp(pattern_minimizers) {
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

    // l is now the lower bound (first index where the suffix is >= pattern)
    l
}


// process_query to handle minimizer-space search, two binary searches, validation, and false positive rate calculation
fn process_query(
    data: &MinimizerStringData,
    original_reference: &str, // Need original reference for validation
    original_query_sequence: String, // The original query DNA sequence
    query_name: &str,
    output_file: &mut File,
) -> () {
    // Transform the original query sequence into its minimizer sequence
    let query_minimizers = compute_minimizers_for_query(
        &original_query_sequence,
        data.minimizer_k,
        data.window_w,
    );

    // If the query minimizer sequence is empty, it cannot match anything.
    if query_minimizers.is_empty() {
         let output_string = format!("{}\t0\t0.0", query_name); // Report 0 occurrences and 0.0 false positive rate
         writeln!(output_file, "{}", output_string).expect("failed to write to the output file!");
         return;
    }

    // --- Find the range [l, r) in the minimizer_sa using two binary searches ---

    // 1. Find the lower bound (l): First suffix >= query_minimizers
    let l = bin_search(
        0, // Search the entire minimizer_sa
        data.minimizer_sa.len(),
        &data.minimizer_sequence,
        &data.minimizer_sa,
        &query_minimizers,
    );

    // 2. Find the upper bound (r): First suffix > query_minimizers
    // Create a pattern that is lexicographically just after query_minimizers.
    // We do this by appending a "terminal minimizer" to the query minimizer sequence.
    let mut upper_bound_pattern = query_minimizers.clone();
    let terminal_minimizer = create_terminal_minimizer_pattern(data.minimizer_k);
    upper_bound_pattern.extend(terminal_minimizer); // Append the terminal minimizer pattern

    let r = bin_search(
        l, // Can start searching from l for efficiency
        data.minimizer_sa.len(),
        &data.minimizer_sequence,
        &data.minimizer_sa,
        &upper_bound_pattern,
    );

    // Calculate the total number of potential matches found in minimizer space
    let total_potential_matches = r.saturating_sub(l);

    // --- Validate potential matches and report original genome positions ---

    let mut true_match_positions: Vec<usize> = Vec::new();

    // Iterate through the range [l, r) in the minimizer_sa
    for i in l..r {
        // Get the corresponding start position in the original genome
        let potential_genome_pos_in_minimizers = data.minimizer_sa[i];

        // Ensure the minimizer sequence match is long enough to correspond
        // to a full query minimizer sequence. A match in the SA means the suffix
        // starts with the query minimizer sequence. We need to check if the
        // reference minimizer sequence at this position is at least as long
        // as the query minimizer sequence.
        if potential_genome_pos_in_minimizers + query_minimizers.len() <= data.minimizer_sequence.len() {

            // Get the start position in the original genome corresponding to the
            // beginning of this minimizer sequence match.
            let original_genome_start_pos = data.minimizer_to_genome_pos[potential_genome_pos_in_minimizers];

            // Perform validation: Check if the original reference sequence
            // at this position matches the original query sequence.
            // Ensure we don't go out of bounds in the original reference.
            if original_genome_start_pos + original_query_sequence.len() <= original_reference.len() {
                let reference_slice = &original_reference[original_genome_start_pos..(original_genome_start_pos + original_query_sequence.len())];

                if reference_slice == original_query_sequence {
                    // This is a true positive match in the original genome
                    true_match_positions.push(original_genome_start_pos);
                }
            }
        }
    }

    // Calculate the number of false positives
    let num_true_matches = true_match_positions.len();
    let num_false_positives = total_potential_matches.saturating_sub(num_true_matches);

    // Calculate the false positive rate
    let false_positive_rate = if total_potential_matches > 0 {
        num_false_positives as f64 / total_potential_matches as f64
    } else {
        0.0
    };


    // Output the results
    let mut output_string = String::new();
    output_string.push_str(query_name);
    output_string.push_str("\t");
    output_string.push_str(&num_true_matches.to_string()); // Report the count of true matches
    output_string.push_str("\t");
    output_string.push_str(&format!("{:.4}", false_positive_rate)); // Report false positive rate, formatted to 4 decimal places

    for pos in true_match_positions {
        output_string.push_str("\t");
        output_string.push_str(&pos.to_string()); // Report the original genome positions
    }

    writeln!(output_file, "{}", output_string).expect("failed to write to the output file!");
}

fn querysa(index: &str, queries: &str, output: &str) -> () {
    // get the serialized reference data
    let data = get_data(index);

    // We need the original reference sequence for validation.
    let original_reference = &data.reference;

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
                    process_query(
                        &data,
                        &original_reference, // Pass the original reference for validation
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
        &original_reference, // Pass the original reference for validation
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
