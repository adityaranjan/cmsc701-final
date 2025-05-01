use bincode;
use std::env;
use std::fs::File;
use std::io::{Read, Write};
use minimizer_sa::shared::{lcp, StringData};

fn get_data(index: &str) -> StringData {
    let mut file = File::open(index).expect("failed to open the file!");
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)
        .expect("failed to read the file!");

    bincode::deserialize(&buffer).expect("failed to deserialize the data!")
}

fn inspectsa(index: &str, sample_rate: usize, output: &str) -> () {
    let StringData {
        reference,
        suffix_array,
        ..
    } = get_data(index);

    let mut mean_lcp: f64 = 0.0;
    let mut max_lcp: f64 = 0.0;
    let mut lcp_list: Vec<f64> = Vec::with_capacity(suffix_array.len() - 1);
    let mut samples: Vec<usize> =
        Vec::with_capacity((suffix_array.len() as f64 / sample_rate as f64).ceil() as usize);

    for i in 0..(suffix_array.len() - 1) {
        let curr_lcp = lcp(
            &reference[suffix_array[i]..],
            &reference[suffix_array[i + 1]..],
        );

        mean_lcp += curr_lcp;

        if curr_lcp > max_lcp {
            max_lcp = curr_lcp;
        }

        lcp_list.push(curr_lcp);

        if i % sample_rate == 0 {
            samples.push(suffix_array[i]);
        }
    }

    if (suffix_array.len() - 1) % sample_rate == 0 {
        samples.push(suffix_array[suffix_array.len() - 1]);
    }

    mean_lcp /= lcp_list.len() as f64;

    lcp_list.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_lcp: f64 = {
        if lcp_list.len() % 2 == 0 {
            (lcp_list[lcp_list.len() / 2] as f64 + lcp_list[(lcp_list.len() / 2) - 1] as f64) / 2.0
        } else {
            lcp_list[lcp_list.len() / 2] as f64
        }
    };
    lcp_list.clear();

    let mut file = File::create(output).expect("failed to create the file!");

    writeln!(file, "{}", mean_lcp).expect("failed to write mean lcp length!");
    writeln!(file, "{}", median_lcp).expect("failed to write median lcp length!");
    writeln!(file, "{}", max_lcp).expect("failed to write maximum lcp length!");
    writeln!(
        file,
        "{}",
        samples
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join("\t")
    )
    .expect("failed to write samples");
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let index = &args[1];

    let sample_rate: usize = args[2].parse().expect("sample_rate must be an integer!");

    let output = &args[3];

    inspectsa(index, sample_rate, output);
}
