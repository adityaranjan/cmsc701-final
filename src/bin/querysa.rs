use minimizer_sa::shared::{compare_minimizer_sequences, compute_minimizers, MinimizerStringData};
use std::cmp::Ordering;
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};

fn get_data(index: &str) -> MinimizerStringData {
    let mut file = File::open(index).expect("failed to open the index file!");
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)
        .expect("failed to read the index file!");

    bincode::deserialize(&buffer).expect("failed to deserialize the index data!")
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
                if m == r - 1 {
                    return r;
                }

                l = m;
            }
            _ => {
                // Suffix is greater than pattern, need to search in the left half
                if m == l + 1 {
                    return m;
                }

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
    partial_check_ct: usize,
    delta_check_ct: usize,
) -> () {
    // Transform the original query sequence into its minimizer sequence (indices into query)
    let mut query_minimizer_indices = compute_minimizers(
        &original_query_sequence,
        data.minimizer_k,
        data.window_w,
        data.minimizer_type,
    );

    // to find lower bound, add a "#" k-mer
    original_query_sequence.push_str(&"#".repeat(data.minimizer_k));
    query_minimizer_indices.push(original_query_sequence.len() - data.minimizer_k);

    // Find the lower bound: First suffix >= query_minimizer_indices
    let l: usize = bin_search(
        0, // Search the entire minimizer_sa
        data.minimizer_sa.len(),
        &data.minimizer_sequence,
        &data.minimizer_sa,
        &query_minimizer_indices,
        &data.reference,
        &original_query_sequence,
        data.minimizer_k,
    );

    // reverse modification made to original query sequence and query minimizer indices
    original_query_sequence.truncate(original_query_sequence.len() - &data.minimizer_k);
    query_minimizer_indices.pop();

    // to find upper bound, add a "}" k-mer
    original_query_sequence.push_str(&"}".repeat(data.minimizer_k));
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

    // reverse modification made to original query sequence and query minimizer indices
    original_query_sequence.truncate(original_query_sequence.len() - &data.minimizer_k);
    query_minimizer_indices.pop();

    // potential matches found in minimizer space

    let potential_match_positions: Vec<(usize, usize)> = (l..r)
        .map(|i| {
            (
                data.minimizer_sequence[data.minimizer_sa[i]],
                data.minimizer_sa[i],
            )
        })
        .collect();

    let mut query_delta = Vec::new();
    for i in 0..(query_minimizer_indices.len() - 1) {
        query_delta.push(query_minimizer_indices[i + 1] - query_minimizer_indices[i]);
    }

    // Output the results
    let mut potential_match_ct = 0;
    let mut output_string = String::new();

    for (pos, minimizer_seq_idx) in potential_match_positions {
        if (pos < query_minimizer_indices[0])
            || (pos + (original_query_sequence.len() - query_minimizer_indices[0] - 1)
                >= data.reference.len() - data.minimizer_k)
        {
            // match isn't valid since the query goes out of bounds
            // even though the minimizers align with the reference
            // note that second condition subtracts minimizer_k to account
            // for the fact that we added a terminal "$" k-mer to the original reference
            continue;
        }

        let mut skip_iter = false;

        // check if deltas of minimizer positions match
        for i in 0..std::cmp::min(delta_check_ct, query_delta.len()) {
            if query_delta[i]
                != (data.minimizer_sequence[minimizer_seq_idx + i + 1]
                    - data.minimizer_sequence[minimizer_seq_idx + i])
            {
                skip_iter = true;
                break;
            }
        }

        if skip_iter {
            // deltas don't match
            continue;
        }

        // check if the minimizer sequence in the reference matches the query
        // for some small number of positions
        let mut remove_indices = HashSet::new();

        for curr_idx in query_minimizer_indices.iter() {
            for i in 0..data.minimizer_k {
                remove_indices.insert(curr_idx + i);
            }
        }

        let mut i = 0;
        let mut relevant_check_ct = 0;

        while relevant_check_ct < partial_check_ct && i < original_query_sequence.len() {
            i += 1;

            if remove_indices.contains(&(i - 1)) {
                continue;
            }

            relevant_check_ct += 1;

            let reference_pos = (i - 1) + pos - query_minimizer_indices[0];

            if &original_query_sequence[(i - 1)..i]
                != &data.reference[reference_pos..(reference_pos + 1)]
            {
                skip_iter = true;
                break;
            }
        }

        if skip_iter {
            // characters don't match
            continue;
        }

        potential_match_ct += 1;

        output_string.push_str("\t");
        output_string.push_str(&((pos - query_minimizer_indices[0]).to_string()));
        // original genome position
    }

    output_string.insert_str(0, &potential_match_ct.to_string());
    output_string.insert_str(0, "\t");
    output_string.insert_str(0, query_name);

    writeln!(output_file, "{}", output_string).expect("failed to write to the output file!");
}

fn querysa(
    index: &str,
    queries: &str,
    output: &str,
    partial_check_ct: usize,
    delta_check_ct: usize,
) -> () {
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
                        partial_check_ct,
                        delta_check_ct,
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
        partial_check_ct,
        delta_check_ct,
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 6 {
        eprintln!(
            "Usage: {} <index_path> <queries_path> <output_path> <partial_check_ct> <delta_check_ct>",
            args[0]
        );
        return;
    }

    let index = &args[1];
    let queries = &args[2];
    let output = &args[3];
    let partial_check_ct: usize = args[4].parse().unwrap_or(0);
    let delta_check_ct: usize = args[5].parse().unwrap_or(0);

    querysa(index, queries, output, partial_check_ct, delta_check_ct);
}
