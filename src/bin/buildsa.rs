use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use minimizer_sa::shared::StringData;

fn get_reference(reference_path: &str) -> String {
    let file = File::open(reference_path).expect("file couldn't be opened!");
    let reader = BufReader::new(file);

    let mut result = reader
        .lines()
        .skip(1)
        .map(|line| line.expect("failed to parse line!"))
        .collect::<String>();
    result.push('$');

    result
}

fn buildsa(reference_path: &str, k: usize, output: &str) -> () {
    let reference = get_reference(reference_path);

    let mut suffix_array: Vec<usize> = (0..reference.len()).collect();

    suffix_array.sort_by(|idx1, idx2| (&reference[*idx1..]).cmp(&reference[*idx2..]));

    let mut min_idx = 0;
    let mut max_idx = 0;
    let mut curr_prefix = "";
    let mut curr_suffix_prefix;
    let mut prefix_lookup = HashMap::<String, (usize, usize)>::new();

    for i in 0..reference.len() {
        if suffix_array[i] + k <= reference.len() {
            curr_suffix_prefix = &reference[(suffix_array[i])..(suffix_array[i] + k)];
            if curr_suffix_prefix != curr_prefix {
                if i > 0 {
                    prefix_lookup.insert(curr_prefix.to_string(), (min_idx, max_idx));
                }

                curr_prefix = curr_suffix_prefix;

                min_idx = i;
                max_idx = i;
            } else {
                max_idx += 1;
            }
        }
    }

    prefix_lookup.insert(curr_prefix.to_string(), (min_idx, max_idx));

    let string_data = StringData {
        reference,
        suffix_array,
        prefix_lookup,
        k,
    };

    let mut output_file = File::create(output).expect("output file creation failed!");
    let serialized_data = bincode::serialize(&string_data).expect("failed to serialize!");
    output_file
        .write_all(&serialized_data)
        .expect("writing the data to the output file failed!");
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let reference_path = &args[1];

    let k: usize = args[2].parse().expect("k must be an integer!");

    let output = &args[3];

    buildsa(reference_path, k, output);
}
