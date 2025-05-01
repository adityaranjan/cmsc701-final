use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use minimizer_sa::shared::{lcp, StringData};

fn get_data(index: &str) -> StringData {
    let mut file = File::open(index).expect("failed to open the file!");
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)
        .expect("failed to read the file!");

    bincode::deserialize(&buffer).expect("failed to deserialize the data!")
}

fn querysa(index: &str, queries: &str, query_mode: &str, output: &str) -> () {
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
                    process_query(
                        &data,
                        curr_sequence,
                        &curr_query,
                        query_mode,
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

    curr_sequence = curr_sequence_vec.join("");
    curr_sequence_vec.clear();
    process_query(
        &data,
        curr_sequence,
        &curr_query,
        query_mode,
        &mut output_file,
    );
}

fn lcp_bin_search(
    left: usize,
    right: usize,
    reference: &str,
    suffix_array: &Vec<usize>,
    pattern: &str,
    cmp_cnt: &mut usize,
    lcp_left: usize,
    lcp_right: usize,
) -> usize {
    let mut l = left;
    let mut r = right;

    let mut lcp_l = lcp_left;
    let mut lcp_r = lcp_right;

    while l < r {
        let m = (l + r) / 2;

        let min_lcp = lcp_l.min(lcp_r);
        let lcp_m = lcp(
            &reference[(suffix_array[m] + min_lcp)..],
            &pattern[min_lcp..],
        ) as usize;

        *cmp_cnt += lcp_m + 1;

        if &pattern[(min_lcp + lcp_m)..] < &reference[(suffix_array[m] + min_lcp + lcp_m)..] {
            if m == l + 1 {
                return m;
            }

            r = m;
            lcp_r = min_lcp + lcp_m;
        } else {
            if m == r - 1 {
                return r;
            }

            l = m;
            lcp_l = min_lcp + lcp_m;
        }
    }

    return l;
}

fn bin_search(
    left: usize,
    right: usize,
    reference: &str,
    suffix_array: &Vec<usize>,
    pattern: &str,
    cmp_cnt: &mut usize,
) -> usize {
    let mut l = left;
    let mut r = right;

    while l < r {
        let m = (l + r) / 2;

        *cmp_cnt += lcp(&reference[suffix_array[m]..], pattern) as usize + 1;

        if &reference[suffix_array[m]..] < pattern {
            if m == r - 1 {
                return r;
            }

            l = m;
        } else {
            if m == l + 1 {
                return m;
            }

            r = m;
        }
    }

    return l;
}

fn process_query(
    data: &StringData,
    curr_sequence: String,
    query_name: &str,
    query_mode: &str,
    output_file: &mut File,
) -> () {
    let mut char_cmp_lb: usize = 0;
    let mut char_cmp_ub: usize = 0;

    let mut pattern = curr_sequence.to_string();

    let mut left_idx = 0;
    let mut right_idx = data.suffix_array.len();

    let k_prefix = &curr_sequence[0..(curr_sequence.len().min((&data).k))];

    let mut pref_lb: usize = 0;
    let mut pref_ub: usize = 0;
    let mut exists_in_lookup = false;
    if query_mode == "prefaccel" {
        if let Some((lb, ub)) = data.prefix_lookup.get(k_prefix) {
            pref_lb = *lb;
            pref_ub = *ub + 1;

            left_idx = pref_lb;
            right_idx = pref_ub;

            exists_in_lookup = true;
        }
    }

    pattern.push('#');

    if exists_in_lookup {
        left_idx = (left_idx - 1).min(0);
    }

    let l: usize;
    let r: usize;

    if query_mode == "simpaccel" {
        l = lcp_bin_search(
            left_idx,
            right_idx,
            &data.reference,
            &data.suffix_array,
            &pattern,
            &mut char_cmp_lb,
            0,
            0,
        );
    } else {
        l = bin_search(
            left_idx,
            right_idx,
            &data.reference,
            &data.suffix_array,
            &pattern,
            &mut char_cmp_lb,
        );
    }

    pattern.pop();
    pattern.push('{');

    if query_mode == "simpaccel" {
        r = lcp_bin_search(
            left_idx,
            right_idx,
            &data.reference,
            &data.suffix_array,
            &pattern,
            &mut char_cmp_ub,
            0,
            0,
        );
    } else {
        r = bin_search(
            left_idx,
            right_idx,
            &data.reference,
            &data.suffix_array,
            &pattern,
            &mut char_cmp_ub,
        );
    }

    let mut output_string = String::new();

    output_string.push_str(query_name);
    output_string.push_str("\t");

    if query_mode == "prefaccel" {
        output_string.push_str(&pref_lb.to_string());
        output_string.push_str("\t");

        output_string.push_str(&pref_ub.to_string());
        output_string.push_str("\t");
    } else {
        output_string.push_str(&char_cmp_lb.to_string());
        output_string.push_str("\t");

        output_string.push_str(&char_cmp_ub.to_string());
        output_string.push_str("\t");
    }

    if l == r {
        output_string.push('0');
    } else {
        output_string.push_str(&((r - l).to_string()));

        for i in l..r {
            output_string.push_str("\t");
            output_string.push_str(&data.suffix_array[i].to_string());
        }
    }

    writeln!(output_file, "{}", output_string).expect("failed to write to the output file!");
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let index = &args[1];
    let queries = &args[2];
    let query_mode = &args[3];
    let output = &args[4];

    querysa(index, queries, query_mode, output);
}
