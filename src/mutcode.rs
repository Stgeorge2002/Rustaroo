//! MutCode edit lists for every member sequence against its family reference.
//!
//! An edit is MutCode's `[start, delete_length, inserted_sequence]`: `start`
//! is a 0-based offset in the original reference, and applying the list
//! rebuilds the member sequence. The gene id is `family/sequence_id` so
//! several copies of one family stay distinct. `call` is `annotated` or
//! `refound`, plus `reference` on a sequence that was selected as the
//! reference.

use anyhow::Result;
use rayon::prelude::*;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::io::Write;

pub const HEADER: &str = "\
sample\tgene\treference_length\treference_sha256\t\
target_length\ttarget_sha256\tcodes\tcall";

pub struct Copy<'a> {
    pub sample: &'a str,
    pub id: &'a str,
    pub dna: &'a str,
    pub refound: bool,
    pub is_reference: bool,
}

#[derive(Debug, Serialize)]
struct Edit(usize, usize, String);

pub fn write_family(
    out: &mut impl Write,
    family: &str,
    reference: &str,
    copies: &[Copy<'_>],
) -> Result<()> {
    let reference_bytes = reference.as_bytes();
    let reference_sha = sha256(reference_bytes);
    let encoded: Vec<_> = copies
        .par_iter()
        .map(|copy| encode(reference_bytes, &reference_sha, copy))
        .collect();
    for (copy, (codes, target_sha)) in copies.iter().zip(encoded) {
        let mut call = if copy.refound { "refound" } else { "annotated" }.to_string();
        if copy.is_reference {
            call.push_str(",reference");
        }
        writeln!(
            out,
            "{}\t{family}/{}\t{}\t{reference_sha}\t{}\t{target_sha}\t{codes}\t{call}",
            copy.sample,
            copy.id,
            reference_bytes.len(),
            copy.dna.len(),
        )?;
    }
    Ok(())
}

fn encode(reference: &[u8], reference_sha: &str, copy: &Copy<'_>) -> (String, String) {
    let target = copy.dna.as_bytes();
    let edits = diff(reference, target);
    let codes = serde_json::to_string(&edits).expect("an edit list is valid JSON");
    debug_assert_eq!(apply(reference, &edits).as_deref(), Some(target));
    let target_sha = if edits.is_empty() {
        reference_sha.to_string()
    } else {
        sha256(target)
    };
    (codes, target_sha)
}

fn sha256(sequence: &[u8]) -> String {
    format!("{:x}", Sha256::digest(sequence))
}

/// Shortest edit script, left-aligned. Substitutions are preferred over an
/// indel pair of the same cost, then deletions over insertions.
fn diff(reference: &[u8], target: &[u8]) -> Vec<Edit> {
    if reference == target {
        return Vec::new();
    }
    let limit = reference.len().max(target.len());
    let mut width = reference.len().abs_diff(target.len()).max(1);
    loop {
        if let Some(edits) = banded(reference, target, width) {
            return edits;
        }
        if width >= limit {
            // A single trimmed replacement always rebuilds the target. The band
            // covers every edit of cost <= limit, so this is only a backstop.
            return banded(reference, target, limit)
                .unwrap_or_else(|| trim(0, reference, target).into_iter().collect());
        }
        width = (width.saturating_mul(2)).min(limit);
    }
}

const MATCH: u8 = 0;
const SUBSTITUTION: u8 = 1;
const DELETION: u8 = 2;
const INSERTION: u8 = 3;

fn banded(reference: &[u8], target: &[u8], width: usize) -> Option<Vec<Edit>> {
    let rows = reference.len();
    let columns = target.len();
    let bounds: Vec<(usize, usize)> = (0..=rows)
        .map(|row| row_bounds(row, rows, columns, width))
        .collect();
    let (end_lo, end_hi) = bounds[rows];
    if columns < end_lo || columns > end_hi {
        return None;
    }

    let mut previous = vec![u32::MAX / 4; columns + 1];
    let mut current = vec![u32::MAX / 4; columns + 1];
    let (_, origin_hi) = bounds[0];
    previous[0] = 0;
    for column in 1..=origin_hi {
        previous[column] = column as u32;
    }

    let mut trace = Vec::with_capacity(rows);
    for row in 1..=rows {
        let (lo, hi) = bounds[row];
        let (previous_lo, previous_hi) = bounds[row - 1];
        current[lo..=hi].fill(u32::MAX / 4);
        let mut directions = vec![DELETION; hi - lo + 1];
        for column in lo..=hi {
            let mut best = u32::MAX / 4;
            let mut direction = DELETION;
            if column > 0 && (column - 1) >= previous_lo && (column - 1) <= previous_hi {
                let substitution =
                    previous[column - 1] + u32::from(reference[row - 1] != target[column - 1]);
                if substitution < best {
                    best = substitution;
                    direction = if reference[row - 1] == target[column - 1] {
                        MATCH
                    } else {
                        SUBSTITUTION
                    };
                }
            }
            if column >= previous_lo && column <= previous_hi {
                let deletion = previous[column].saturating_add(1);
                if deletion < best {
                    best = deletion;
                    direction = DELETION;
                }
            }
            if column > 0 && (column - 1) >= lo && (column - 1) <= hi {
                let insertion = current[column - 1].saturating_add(1);
                if insertion < best {
                    best = insertion;
                    direction = INSERTION;
                }
            }
            current[column] = best;
            directions[column - lo] = direction;
        }
        trace.push(directions);
        std::mem::swap(&mut previous, &mut current);
    }

    let cost = previous[columns];
    if cost as usize > width {
        return None;
    }

    let mut operations = Vec::with_capacity(rows + columns);
    let mut row = rows;
    let mut column = columns;
    while row > 0 || column > 0 {
        if row == 0 {
            operations.push(INSERTION);
            column -= 1;
            continue;
        }
        if column == 0 {
            operations.push(DELETION);
            row -= 1;
            continue;
        }
        let (lo, hi) = bounds[row];
        if column < lo || column > hi {
            return None;
        }
        match trace[row - 1][column - lo] {
            DELETION => {
                operations.push(DELETION);
                row -= 1;
            }
            INSERTION => {
                operations.push(INSERTION);
                column -= 1;
            }
            kind => {
                operations.push(kind);
                row -= 1;
                column -= 1;
            }
        }
    }
    operations.reverse();
    Some(edits_of(reference, target, &operations))
}

fn row_bounds(row: usize, rows: usize, columns: usize, width: usize) -> (usize, usize) {
    let row = row as isize;
    let rows = rows as isize;
    let columns = columns as isize;
    let width = width as isize;
    let difference = columns - rows;
    let low = (row + difference.min(0) - width).max(0) as usize;
    let high = ((row + difference.max(0) + width).max(0) as usize).min(columns as usize);
    (low, high)
}

fn edits_of(reference: &[u8], target: &[u8], operations: &[u8]) -> Vec<Edit> {
    let mut edits = Vec::new();
    let mut reference_index = 0;
    let mut target_index = 0;
    let mut cursor = 0;
    while cursor < operations.len() {
        if operations[cursor] == MATCH {
            reference_index += 1;
            target_index += 1;
            cursor += 1;
            continue;
        }
        let start = reference_index;
        let inserted_at = target_index;
        while cursor < operations.len() && operations[cursor] != MATCH {
            match operations[cursor] {
                DELETION => reference_index += 1,
                INSERTION => target_index += 1,
                _ => {
                    reference_index += 1;
                    target_index += 1;
                }
            }
            cursor += 1;
        }
        if let Some(edit) = trim(
            start,
            &reference[start..reference_index],
            &target[inserted_at..target_index],
        ) {
            edits.push(edit);
        }
    }
    edits
}

/// MutCode's prefix/suffix trim of one already-placed replacement.
fn trim(start: usize, deleted: &[u8], inserted: &[u8]) -> Option<Edit> {
    let prefix = deleted
        .iter()
        .zip(inserted)
        .take_while(|(reference, alternate)| reference == alternate)
        .count();
    let suffix = deleted[prefix..]
        .iter()
        .rev()
        .zip(inserted[prefix..].iter().rev())
        .take_while(|(reference, alternate)| reference == alternate)
        .count();
    let delete_length = deleted.len() - prefix - suffix;
    let insert_end = inserted.len() - suffix;
    let insert = &inserted[prefix..insert_end];
    if delete_length == 0 && insert.is_empty() {
        None
    } else {
        Some(Edit(
            start + prefix,
            delete_length,
            String::from_utf8(insert.to_vec()).expect("member DNA is ASCII"),
        ))
    }
}

fn apply(reference: &[u8], edits: &[Edit]) -> Option<Vec<u8>> {
    let mut sequence = Vec::new();
    let mut cursor = 0;
    for Edit(start, deleted, inserted) in edits {
        if *start < cursor || start + deleted > reference.len() {
            return None;
        }
        sequence.extend_from_slice(&reference[cursor..*start]);
        sequence.extend_from_slice(inserted.as_bytes());
        cursor = start + deleted;
    }
    sequence.extend_from_slice(&reference[cursor..]);
    Some(sequence)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_match_mutcode_edits_and_round_trip() {
        let reference = "ATCGATGCATAGGCTTAC";
        let mut target = reference.to_string();
        target.replace_range(4..5, "T");
        target.insert_str(9, "CC");
        target.replace_range(15..16, "");
        let edits = diff(reference.as_bytes(), target.as_bytes());
        assert_eq!(
            serde_json::to_string(&edits).unwrap(),
            r#"[[4,1,"T"],[9,0,"CC"],[13,1,""]]"#
        );
        assert_eq!(
            apply(reference.as_bytes(), &edits).unwrap(),
            target.as_bytes()
        );
    }

    #[test]
    fn indels_are_left_aligned() {
        assert_eq!(
            serde_json::to_string(&diff(b"AAAAA", b"AAAA")).unwrap(),
            r#"[[0,1,""]]"#
        );
        assert_eq!(
            serde_json::to_string(&diff(b"AAAA", b"AAAAA")).unwrap(),
            r#"[[0,0,"A"]]"#
        );
        assert_eq!(
            serde_json::to_string(&diff(b"AAACCCGGGTTT", b"AATCCCAGGTTT")).unwrap(),
            r#"[[2,1,"T"],[6,1,"A"]]"#
        );
    }

    #[test]
    fn identical_sequence_has_an_empty_list() {
        let sequence = b"ATGAAA";
        assert!(diff(sequence, sequence).is_empty());
        assert_eq!(
            serde_json::to_string(&diff(sequence, sequence)).unwrap(),
            "[]"
        );
    }

    #[test]
    fn call_column_marks_the_selected_reference() {
        let mut output = Vec::new();
        write_family(
            &mut output,
            "tetM",
            "ATGC",
            &[
                Copy {
                    sample: "sample_12",
                    id: "tetM_0001",
                    dna: "ATGC",
                    refound: false,
                    is_reference: true,
                },
                Copy {
                    sample: "sample_12",
                    id: "tetM_0002",
                    dna: "ATGG",
                    refound: false,
                    is_reference: false,
                },
                Copy {
                    sample: "sample_03",
                    id: "tetM_refound_1",
                    dna: "ATG",
                    refound: true,
                    is_reference: false,
                },
            ],
        )
        .unwrap();
        let text = String::from_utf8(output).unwrap();
        let mut lines = text.lines();
        assert!(lines
            .next()
            .unwrap()
            .starts_with("sample_12\ttetM/tetM_0001\t"));
        assert!(lines.next().unwrap().ends_with("\tannotated"));
        assert!(text
            .lines()
            .next()
            .unwrap()
            .ends_with("\tannotated,reference"));
        assert!(text.lines().nth(1).unwrap().contains("\t[[3,1,\"G\"]]\t"));
        assert!(text.lines().nth(2).unwrap().ends_with("\trefound"));
    }

    /// Same tie order as `banded`, on the whole matrix, so a band that drops
    /// an optimal path cannot hide.
    fn full_diff(reference: &[u8], target: &[u8]) -> Vec<Edit> {
        let rows = reference.len();
        let columns = target.len();
        let mut cost = vec![vec![0u32; columns + 1]; rows + 1];
        let mut trace = vec![vec![DELETION; columns + 1]; rows + 1];
        for column in 1..=columns {
            cost[0][column] = column as u32;
            trace[0][column] = INSERTION;
        }
        for row in 1..=rows {
            cost[row][0] = row as u32;
            for column in 1..=columns {
                let substitution =
                    cost[row - 1][column - 1] + u32::from(reference[row - 1] != target[column - 1]);
                let deletion = cost[row - 1][column] + 1;
                let insertion = cost[row][column - 1] + 1;
                let (best, direction) = if substitution <= deletion && substitution <= insertion {
                    (
                        substitution,
                        if reference[row - 1] == target[column - 1] {
                            MATCH
                        } else {
                            SUBSTITUTION
                        },
                    )
                } else if deletion <= insertion {
                    (deletion, DELETION)
                } else {
                    (insertion, INSERTION)
                };
                cost[row][column] = best;
                trace[row][column] = direction;
            }
        }
        let mut operations = Vec::new();
        let mut row = rows;
        let mut column = columns;
        while row > 0 || column > 0 {
            match trace[row][column] {
                DELETION => {
                    operations.push(DELETION);
                    row -= 1;
                }
                INSERTION => {
                    operations.push(INSERTION);
                    column -= 1;
                }
                kind => {
                    operations.push(kind);
                    row -= 1;
                    column -= 1;
                }
            }
        }
        operations.reverse();
        edits_of(reference, target, &operations)
    }

    #[test]
    fn band_matches_the_full_matrix_for_short_sequences() {
        let mut sequences = vec![Vec::new()];
        for length in 1..=4 {
            let mut grown = Vec::new();
            for sequence in &sequences {
                if sequence.len() + 1 != length {
                    continue;
                }
                for base in [b'A', b'C', b'G', b'T'] {
                    let mut next = sequence.clone();
                    next.push(base);
                    grown.push(next);
                }
            }
            sequences.extend(grown);
        }
        for reference in &sequences {
            for target in &sequences {
                let band = diff(reference, target);
                let full = full_diff(reference, target);
                assert_eq!(
                    serde_json::to_string(&band).unwrap(),
                    serde_json::to_string(&full).unwrap(),
                    "{reference:?} -> {target:?}"
                );
                assert_eq!(apply(reference, &band).unwrap(), target.as_slice());
            }
        }
    }

    /// Gene-length pairs, checked with MutCode's own rules: sorted,
    /// non-overlapping edits, IUPAC inserts, and byte-for-byte apply.
    #[test]
    fn gene_length_edits_follow_mutcode_rules() {
        let mut state = 0x1234_5678_u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1);
            state
        };
        let bases = [b'A', b'C', b'G', b'T'];
        for _ in 0..30 {
            let length = 80 + (next() as usize % 240);
            let reference: Vec<u8> = (0..length).map(|_| bases[(next() % 4) as usize]).collect();
            let mut target = reference.clone();
            let changes = 1 + (next() as usize % 12);
            for _ in 0..changes {
                let at = (next() as usize) % target.len().max(1);
                match next() % 3 {
                    0 => target[at] = bases[(next() % 4) as usize],
                    1 => {
                        target.insert(at, bases[(next() % 4) as usize]);
                    }
                    _ => {
                        target.remove(at);
                    }
                }
                if target.is_empty() {
                    target.push(b'A');
                }
            }
            let edits = diff(&reference, &target);
            let mut cursor = 0;
            let mut previous = None;
            for edit in &edits {
                let Edit(start, deleted, inserted) = edit;
                let end = start + deleted;
                assert!(*start >= cursor && end <= reference.len(), "{edit:?}");
                assert!(previous.is_none_or(|previous: usize| previous < *start));
                assert!(!(*deleted == 0 && inserted.is_empty()));
                assert!(inserted.bytes().all(|base| matches!(
                    base,
                    b'A' | b'C' | b'G' | b'T' | b'N'
                )));
                cursor = end;
                previous = Some(*start);
            }
            assert_eq!(apply(&reference, &edits).unwrap(), target);
            assert_eq!(
                serde_json::to_string(&edits).unwrap(),
                serde_json::to_string(&full_diff(&reference, &target)).unwrap()
            );
        }
    }
}
