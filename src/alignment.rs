//! Sequence alignment utilities using edlib-style algorithms

use bio::pattern_matching::myers::long::Myers;

#[derive(Default)]
pub struct HwEditScratch {
    previous: Vec<usize>,
    current: Vec<usize>,
}

impl HwEditScratch {
    /// Compute the exact HW edit distance using reusable linear working memory.
    pub fn distance(&mut self, query: &[u8], target: &[u8]) -> usize {
        self.distance_with_limit(query, target, usize::MAX)
            .unwrap_or(usize::MAX)
    }

    /// Return the exact distance when it is within the limit, otherwise `None`.
    pub fn distance_with_limit(
        &mut self,
        query: &[u8],
        target: &[u8],
        max_distance: usize,
    ) -> Option<usize> {
        let m = query.len();
        let n = target.len();
        if m == 0 || n == 0 {
            return None;
        }

        if m <= n && target.windows(m).any(|window| window == query) {
            return Some(0);
        }

        let cap = max_distance.saturating_add(1);
        self.previous.resize(n + 1, 0);
        self.current.resize(n + 1, 0);
        self.previous.fill(0);

        for i in 1..=m {
            self.current[0] = i.min(cap);
            for j in 1..=n {
                let substitution = self.previous[j - 1]
                    .saturating_add(usize::from(query[i - 1] != target[j - 1]))
                    .min(cap);
                self.current[j] = substitution
                    .min(self.previous[j].saturating_add(1).min(cap))
                    .min(self.current[j - 1].saturating_add(1).min(cap));
            }
            std::mem::swap(&mut self.previous, &mut self.current);
            let row_minimum = self.previous.iter().copied().min().unwrap_or(cap);
            if row_minimum > max_distance {
                return None;
            }
        }

        let distance = self.previous[1..].iter().copied().min().unwrap_or(cap);
        (distance <= max_distance).then_some(distance)
    }
}

/// Compute the exact HW edit distance using linear working memory.
pub fn hw_edit_distance(query: &[u8], target: &[u8]) -> usize {
    HwEditScratch::default().distance(query, target)
}

/// Edlib-style alignment with HW (prefix-suffix free) mode
pub fn align_sequences(query: &[u8], target: &[u8], max_distance: usize) -> AlignmentResult {
    let m = query.len();
    let n = target.len();
    let no_match = AlignmentResult {
        edit_distance: usize::MAX,
        locations: vec![],
        alignment: None,
    };
    if m == 0 || n == 0 {
        return no_match;
    }

    let matcher = Myers::<u64>::new(query);
    let mut best_dist = usize::MAX;
    let mut best_end = 0;
    for (end, distance) in matcher.find_all_end(target, max_distance.min(m)) {
        if distance < best_dist {
            best_dist = distance;
            best_end = end + 1;
            if distance == 0 {
                break;
            }
        }
    }
    if best_dist == usize::MAX {
        return no_match;
    }
    if best_dist == 0 {
        return AlignmentResult {
            edit_distance: 0,
            locations: vec![(best_end - m, best_end)],
            alignment: None,
        };
    }

    // An optimal path consumes at most m + distance target bases. Propagate its
    // origin in that window, retaining the old substitution/deletion/insertion ties.
    let offset = best_end.saturating_sub(m.saturating_add(best_dist));
    let window = &target[offset..best_end];
    let width = window.len();
    let mut previous = vec![0; width + 1];
    let mut current = vec![0; width + 1];
    let mut previous_start: Vec<usize> = (0..=width).collect();
    let mut current_start = vec![0; width + 1];
    for i in 1..=m {
        current[0] = i;
        current_start[0] = 0;
        for j in 1..=width {
            let substitution = previous[j - 1] + usize::from(query[i - 1] != window[j - 1]);
            let deletion = previous[j] + 1;
            let insertion = current[j - 1] + 1;
            current[j] = substitution.min(deletion).min(insertion);
            current_start[j] = if current[j] == substitution {
                previous_start[j - 1]
            } else if current[j] == deletion {
                previous_start[j]
            } else {
                current_start[j - 1]
            };
        }
        std::mem::swap(&mut previous, &mut current);
        std::mem::swap(&mut previous_start, &mut current_start);
    }
    debug_assert_eq!(previous[width], best_dist);

    AlignmentResult {
        edit_distance: best_dist,
        locations: vec![(offset + previous_start[width], best_end)],
        alignment: None,
    }
}

#[derive(Debug, Clone)]
pub struct AlignmentResult {
    pub edit_distance: usize,
    pub locations: Vec<(usize, usize)>,
    pub alignment: Option<String>,
}

pub fn calculate_identity(edit_distance: usize, alignment_length: usize) -> f64 {
    if alignment_length == 0 {
        return 0.0;
    }
    1.0 - (edit_distance as f64 / alignment_length as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn full_table_alignment(query: &[u8], target: &[u8], limit: usize) -> AlignmentResult {
        let m = query.len();
        let n = target.len();
        let mut result = AlignmentResult {
            edit_distance: usize::MAX,
            locations: Vec::new(),
            alignment: None,
        };
        if m == 0 || n == 0 {
            return result;
        }
        let mut dp = vec![vec![0; n + 1]; m + 1];
        let mut traceback = vec![vec![(0, 0); n + 1]; m + 1];
        for i in 0..=m {
            dp[i][0] = i;
        }
        for i in 1..=m {
            for j in 1..=n {
                let substitution = dp[i - 1][j - 1] + usize::from(query[i - 1] != target[j - 1]);
                let deletion = dp[i - 1][j] + 1;
                let insertion = dp[i][j - 1] + 1;
                dp[i][j] = substitution.min(deletion).min(insertion);
                traceback[i][j] = if dp[i][j] == substitution {
                    (i - 1, j - 1)
                } else if dp[i][j] == deletion {
                    (i - 1, j)
                } else {
                    (i, j - 1)
                };
            }
        }
        let mut end = 0;
        for j in 1..=n {
            if dp[m][j] < result.edit_distance {
                result.edit_distance = dp[m][j];
                end = j;
            }
        }
        if result.edit_distance > limit {
            result.edit_distance = usize::MAX;
            return result;
        }
        let (mut i, mut j) = (m, end);
        while i > 0 {
            (i, j) = traceback[i][j];
        }
        result.locations.push((j, end));
        result
    }

    #[test]
    fn bit_parallel_alignment_preserves_distances_locations_and_ties() {
        let mut sequences = vec![Vec::new()];
        for length in 1..=5 {
            for bits in 0..(1 << length) {
                sequences.push(
                    (0..length)
                        .map(|bit| if bits & (1 << bit) == 0 { b'A' } else { b'C' })
                        .collect(),
                );
            }
        }
        for query in &sequences {
            for target in &sequences {
                for limit in [0, 1, 2, 4, usize::MAX] {
                    let expected = full_table_alignment(query, target, limit);
                    let actual = align_sequences(query, target, limit);
                    assert_eq!(
                        (actual.edit_distance, actual.locations),
                        (expected.edit_distance, expected.locations),
                        "query={query:?}, target={target:?}, limit={limit}"
                    );
                }
            }
        }
    }

    #[test]
    fn bit_parallel_alignment_preserves_multiblock_indels() {
        for length in [63, 64, 65, 127, 128, 129, 257] {
            let query: Vec<_> = (0..length).map(|i| b"ACGT"[i % 4]).collect();
            let mut target = b"TTTTTT".to_vec();
            target.extend(&query);
            target.insert(length / 2, b'A');
            target.remove(length / 3);
            target[length / 4] = b'N';
            target.extend(b"TTTTTT");
            for limit in [0, 1, 2, 3, usize::MAX] {
                let expected = full_table_alignment(&query, &target, limit);
                let actual = align_sequences(&query, &target, limit);
                assert_eq!(
                    (actual.edit_distance, actual.locations),
                    (expected.edit_distance, expected.locations)
                );
            }
        }
    }

    #[test]
    fn bit_parallel_alignment_searches_long_contigs() {
        let query = vec![b'C'; 96];
        let mut target = vec![b'G'; 100_000];
        target.extend(&query);
        target.extend(&query);
        let actual = align_sequences(&query, &target, 0);
        assert_eq!(actual.edit_distance, 0);
        assert_eq!(actual.locations, vec![(100_000, 100_096)]);
    }

    #[test]
    fn test_exact_match() {
        let query = b"ACGT";
        let target = b"ACGT";
        let result = align_sequences(query, target, 10);
        assert_eq!(result.edit_distance, 0);
    }

    #[test]
    fn test_with_mismatches() {
        let query = b"ACGT";
        let target = b"ACTT";
        let result = align_sequences(query, target, 10);
        assert!(result.edit_distance <= 1);
    }

    #[test]
    fn test_linear_space_hw_distance_matches_full_alignment() {
        let sequences = [
            b"A".as_slice(),
            b"C".as_slice(),
            b"AA".as_slice(),
            b"AC".as_slice(),
            b"CA".as_slice(),
            b"CC".as_slice(),
            b"AAA".as_slice(),
            b"AAC".as_slice(),
            b"ACA".as_slice(),
            b"ACC".as_slice(),
        ];

        for query in sequences {
            for target in sequences {
                let mut scratch = HwEditScratch::default();
                let expected = align_sequences(query, target, usize::MAX).edit_distance;
                assert_eq!(scratch.distance(query, target), expected);
                for limit in 0..=3 {
                    let bounded = scratch.distance_with_limit(query, target, limit);
                    let expected = (expected <= limit).then_some(expected);
                    assert_eq!(bounded, expected);
                }
            }
        }
    }
}
