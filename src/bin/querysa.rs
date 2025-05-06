use std::env;
use std::fs::File;
use std::cmp::Ordering;
use std::io::{BufRead, BufReader, Read, Write};
use minimizer_sa::shared::{MinimizerStringData, compute_minimizers, compare_minimizer_sequences};


fn get_data(index: &str) -> MinimizerStringData {
    let mut file = File::open(index).expect("failed to open the index file!");
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)
        .expect("failed to read the index file!");

    bincode::deserialize(&buffer).expect("failed to deserialize the index data!")
}


fn create_terminal_minimizer_pattern(k: usize) -> String {
    // A single k-mer string of '{' characters, which is lexicographically
    // greater than any k-mer composed of A, C, G, T.
    "{".repeat(k)
}


// Modified bin_search to work on sequences of original positions.
fn bin_search(
    left: usize,
    right: usize,
    minimizer_sequence_indices: &Vec<usize>,
    minimizer_sa: &Vec<usize>,
    pattern_minimizer_indices: &Vec<usize>,
    original_reference: &str,
    original_query_sequence: &str,
    k: usize,
) -> usize {
    let mut l = left;
    let mut r = right;

    while l < r {
        let m = l + (r - l) / 2;

        // Get the slice of original positions for the suffix in the reference minimizer sequence
        let suffix_in_minimizers_indices = &minimizer_sequence_indices[minimizer_sa[m]..];

        // Compare the suffix (sequence of indices from reference) with the pattern (sequence of indices from query)
        match compare_minimizer_sequences(
            None,
            None,
            suffix_in_minimizers_indices,
            pattern_minimizer_indices,
            original_reference,
            original_query_sequence,
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

    l
}


fn process_query(
    data: &MinimizerStringData,
    mut original_query_sequence: String,
    query_name: &str,
    output_file: &mut File,
) -> () {
    // Transform the original query sequence into its minimizer sequence (indices into query)
    let mut query_minimizer_indices = compute_minimizers(
        &original_query_sequence,
        data.minimizer_k,
        data.window_w,
    );

    // Find the lower bound: First suffix >= query_minimizer_indices
    let l = bin_search(
        0, // Search the entire minimizer_sa
        data.minimizer_sa.len(),
        &data.minimizer_sequence,
        &data.minimizer_sa,
        &query_minimizer_indices,
        &data.reference,
        &original_query_sequence,
        data.minimizer_k,
    );

    // add terminal minimizer to query
    original_query_sequence.push_str(&create_terminal_minimizer_pattern(data.minimizer_k));
    query_minimizer_indices.push(original_query_sequence.len() - data.minimizer_k);

    // Find the upper bound: First suffix > query_minimizer_indices
    let r = bin_search(
        l, // Start at l for efficiency
        data.minimizer_sa.len(),
        &data.minimizer_sequence,
        &data.minimizer_sa,
        &query_minimizer_indices,
        &data.reference,
        &original_query_sequence,
        data.minimizer_k,
    );

    // potential matches found in minimizer space

    let potential_match_positions: Vec<usize> = (l..r)
        .map(|i| data.minimizer_sequence[data.minimizer_sa[i]]) 
        .collect();

    // Output the results
    let mut output_string = String::new();
    output_string.push_str(query_name);
    output_string.push_str("\t");
    output_string.push_str(&potential_match_positions.len().to_string()); // count of potential matches

    for pos in potential_match_positions {
        output_string.push_str("\t");
        output_string.push_str(&pos.to_string()); // original genome positions
    }

    writeln!(output_file, "{}", output_string).expect("failed to write to the output file!");
}

fn querysa(index: &str, queries: &str, output: &str) -> () {
    let data = get_data(index);

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
