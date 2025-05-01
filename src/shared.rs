use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Serialize, Deserialize)]
pub struct StringData {
    pub reference: String,
    pub suffix_array: Vec<usize>,
    pub prefix_lookup: HashMap<String, (usize, usize)>,
    pub k: usize,
}

pub fn lcp(s1: &str, s2: &str) -> f64 {
    let mut i: usize = 0;

    while i < s1.len() && i < s2.len() && &s1[i..i + 1] == &s2[i..i + 1] {
        i += 1;
    }

    i as f64
}
