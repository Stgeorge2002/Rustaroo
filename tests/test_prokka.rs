//! Tests for prokka module
//!
//! Port of tests/test_isvalid.py and related tests

use rustaroo::prokka::{reverse_complement, translate};

#[test]
fn test_translate_basic() {
    // Test basic translation
    let dna = "ATGGCCTAG";
    let protein = translate(dna, 11);
    assert_eq!(protein, "MA*");
}

#[test]
fn test_translate_longer_sequence() {
    let dna = "ATGAAATAG";
    let protein = translate(dna, 11);
    assert_eq!(protein, "MK*");
}

#[test]
fn test_reverse_complement() {
    assert_eq!(reverse_complement("ACGT"), "ACGT");
    assert_eq!(reverse_complement("AAAA"), "TTTT");
    assert_eq!(reverse_complement("ATCG"), "CGAT");
    assert_eq!(reverse_complement("GCTA"), "TAGC");
}

#[test]
fn test_reverse_complement_mixed_case() {
    assert_eq!(reverse_complement("AtCg"), "CgAt");
}
