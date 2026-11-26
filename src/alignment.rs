//! Sequence alignment utilities using edlib-style algorithms

use bio::pattern_matching::myers::long::Myers;
#[cfg(test)]
use bio::pattern_matching::myers::MyersBuilder;
use once_cell::sync::OnceCell;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RefindAlignment {
    pub distance: usize,
    pub start: usize,
    pub end: usize,
    pub terminal_insertions: usize,
}

fn refind_equal(a: u8, b: u8) -> bool {
    a == b
        || (b"ACGT".contains(&a) && b"NE".contains(&b))
        || (b"NE".contains(&a) && b"ACGT".contains(&b))
}

pub(crate) struct RefindMatcher {
    query: Vec<u8>,
    length: usize,
    #[cfg(test)]
    matcher: Myers<u64>,
    forward: OnceCell<GlobalPattern>,
    reverse: OnceCell<GlobalPattern>,
}

#[derive(Clone, Copy)]
struct GlobalState {
    positive: u64,
    negative: u64,
    distance: usize,
}

struct GlobalPattern {
    slots: [u16; 256],
    masks: Vec<u64>,
    boundaries: Vec<u64>,
    /// Bit index of each block's boundary row (63 for every full block).
    shifts: Vec<u32>,
    length: usize,
}

impl GlobalPattern {
    fn new(query: impl ExactSizeIterator<Item = u8> + Clone) -> Self {
        fn matching_bytes(base: u8) -> impl Iterator<Item = u8> {
            let equivalents: &[u8] = match base {
                b'A' | b'C' | b'G' | b'T' => b"NE",
                b'N' | b'E' => b"ACGT",
                _ => b"",
            };
            std::iter::once(base).chain(equivalents.iter().copied())
        }

        let length = query.len();
        let block_count = length.div_ceil(64);
        // Slot zero matches nothing. Other slots retain literal matches even
        // for bytes outside the usual DNA alphabet.
        let mut slots = [0; 256];
        let mut slot_count = 1;
        for symbol in query.clone().flat_map(matching_bytes) {
            let slot = &mut slots[usize::from(symbol)];
            if *slot == 0 {
                *slot = slot_count;
                slot_count += 1;
            }
        }
        let mut masks = vec![0; usize::from(slot_count) * block_count];
        for (i, base) in query.enumerate() {
            let bit = 1u64 << (i % 64);
            for symbol in matching_bytes(base) {
                let slot = usize::from(slots[usize::from(symbol)]);
                masks[slot * block_count + i / 64] |= bit;
            }
        }
        let boundaries: Vec<u64> = (0..block_count)
            .map(|block| 1u64 << ((length - block * 64).min(64) - 1))
            .collect();
        let shifts = boundaries.iter().map(|boundary| boundary.trailing_zeros()).collect();
        Self {
            slots,
            masks,
            boundaries,
            shifts,
            length,
        }
    }

    fn reset(&self, states: &mut Vec<GlobalState>) {
        states.clear();
        for block in 0..self.boundaries.len() {
            states.push(GlobalState {
                positive: u64::MAX,
                negative: 0,
                distance: ((block + 1) * 64).min(self.length),
            });
        }
    }

    fn advance(&self, states: &mut [GlobalState], base: u8) -> usize {
        // A global alignment charges for every target prefix: the top row
        // increases by one, unlike the free-prefix HW endpoint search.
        let block_count = self.boundaries.len();
        let start = usize::from(self.slots[usize::from(base)]) * block_count;
        let masks = &self.masks[start..start + block_count];
        let (mut positive_carry, mut negative_carry) = (1, 0);
        for ((&equal, &shift), state) in masks.iter().zip(&self.shifts).zip(states.iter_mut()) {
            (positive_carry, negative_carry) =
                advance_block_bits(state, equal, shift, positive_carry, negative_carry);
        }
        states.last().expect("nonempty query").distance
    }

    /// `advance` restricted to blocks `low..high`, with the top boundary of
    /// block `low` increasing by one per column; returns the distance of the
    /// last processed block.
    fn advance_range(&self, states: &mut [GlobalState], base: u8, low: usize, high: usize) -> usize {
        let block_count = self.boundaries.len();
        let start = usize::from(self.slots[usize::from(base)]) * block_count;
        let masks = &self.masks[start + low..start + high];
        let (mut positive_carry, mut negative_carry) = (1, 0);
        for ((&equal, &shift), state) in masks
            .iter()
            .zip(&self.shifts[low..high])
            .zip(states[low..high].iter_mut())
        {
            (positive_carry, negative_carry) =
                advance_block_bits(state, equal, shift, positive_carry, negative_carry);
        }
        states[high - 1].distance
    }
}

/// The same block step as `advance_refind_block`, with the incoming and
/// outgoing horizontal delta carried as two 0/1 words (positive, negative)
/// and the boundary row given as a bit index. Feeding the bits straight into
/// the shifts and masks avoids the sign tests and conversions of the `i8`
/// carry on the hot path; the arithmetic, and so every result, is identical.
#[inline(always)]
fn advance_block_bits(
    state: &mut GlobalState,
    equal: u64,
    shift: u32,
    positive_carry: u64,
    negative_carry: u64,
) -> (u64, u64) {
    let vertical = equal | state.negative;
    let equal = equal | negative_carry;
    let horizontal =
        ((equal & state.positive).wrapping_add(state.positive) ^ state.positive) | equal;
    let positive = state.negative | !(horizontal | state.positive);
    let negative = state.positive & horizontal;
    // At most one of the two boundary bits is set (positive and negative
    // horizontal deltas are disjoint).
    let positive_out = (positive >> shift) & 1;
    let negative_out = (negative >> shift) & 1;
    state.distance = state
        .distance
        .wrapping_add(positive_out as usize)
        .wrapping_sub(negative_out as usize);
    let positive = (positive << 1) | positive_carry;
    let negative = (negative << 1) | negative_carry;
    state.positive = negative | !(vertical | positive);
    state.negative = positive & vertical;
    (positive_out, negative_out)
}

#[inline]
fn advance_refind_block(state: &mut GlobalState, equal: u64, boundary: u64, carry: i8) -> i8 {
    let vertical = equal | state.negative;
    let equal = equal | u64::from(carry < 0);
    let horizontal =
        ((equal & state.positive).wrapping_add(state.positive) ^ state.positive) | equal;
    let mut positive = state.negative | !(horizontal | state.positive);
    let mut negative = state.positive & horizontal;
    let next_carry = i8::from(positive & boundary != 0) - i8::from(negative & boundary != 0);
    state.distance = state.distance.wrapping_add_signed(next_carry as isize);
    positive = (positive << 1) | u64::from(carry > 0);
    negative = (negative << 1) | u64::from(carry < 0);
    state.positive = negative | !(vertical | positive);
    state.negative = positive & vertical;
    next_carry
}

#[derive(Default)]
pub(crate) struct RefindScratch {
    states: Vec<GlobalState>,
    columns: Vec<GlobalState>,
    previous: Vec<usize>,
    current: Vec<usize>,
    paths: Vec<u8>,
    row_offsets: Vec<usize>,
    row_starts: Vec<usize>,
}

impl RefindScratch {
    fn distance(&self, pattern: &GlobalPattern, row: usize, column: usize) -> usize {
        if row == 0 {
            return column;
        }
        if column == 0 {
            return row;
        }
        let block = (row - 1) / 64;
        let state = self.columns[(column - 1) * pattern.boundaries.len() + block];
        let local_row = (row - 1) % 64 + 1;
        let above = if local_row == 64 {
            0
        } else {
            u64::MAX << local_row
        };
        let used = u64::MAX >> (64 - (pattern.length - block * 64).min(64));
        let mask = above & used;
        state.distance + (state.negative & mask).count_ones() as usize
            - (state.positive & mask).count_ones() as usize
    }
}

impl RefindMatcher {
    pub fn new(query: &[u8]) -> Option<Self> {
        if query.is_empty() {
            return None;
        }
        Some(Self {
            query: query.to_vec(),
            length: query.len(),
            #[cfg(test)]
            matcher: {
                let mut builder = MyersBuilder::new();
                for base in b"ACGT" {
                    builder.ambig(*base, b"NE");
                }
                builder.ambig(b'N', b"ACGT");
                builder.ambig(b'E', b"ACGT");
                builder.build_long_64(query)
            },
            forward: OnceCell::new(),
            reverse: OnceCell::new(),
        })
    }

    #[cfg(test)]
    pub fn align(&self, target: &[u8]) -> Option<RefindAlignment> {
        self.align_with_scratch(target, &mut RefindScratch::default())
    }

    pub(crate) fn align_with_scratch(
        &self,
        target: &[u8],
        scratch: &mut RefindScratch,
    ) -> Option<RefindAlignment> {
        if target.is_empty() {
            return None;
        }
        let (best, first_end, last_end) = self.single_pass_endpoints(target, scratch)?;
        Some(self.recover(target, best, first_end, last_end, scratch))
    }

    /// A cheap first distance limit: zero when the query occurs verbatim in
    /// the target, since that occurrence is a zero-cost global alignment.
    fn first_limit(&self, target: &[u8]) -> Option<usize> {
        // Only worth it when the limit can prune whole 64-row blocks.
        if self.length <= 64 {
            return None;
        }
        memchr::memmem::find(target, &self.query).map(|_| 0)
    }

    fn single_pass_endpoints(
        &self,
        target: &[u8],
        scratch: &mut RefindScratch,
    ) -> Option<(usize, usize, usize)> {
        // Every cell of the last row is at most the query length, so an
        // unlimited scan always finds an endpoint for a nonempty target. A
        // scan limited to `limit` is exact whenever it finds an endpoint: the
        // block cut-off only discards cells whose distance exceeds the bound,
        // and the bound never drops below the optimum. If the limited scan
        // finds nothing, the optimum exceeds `limit` and the full scan runs
        // (the verbatim-occurrence limit of zero always succeeds).
        if let Some(limit) = self.first_limit(target) {
            if let Some(found) = self.single_pass_endpoints_within(target, scratch, limit) {
                return Some(found);
            }
        }
        self.single_pass_endpoints_within(target, scratch, self.length)
    }

    fn single_pass_endpoints_within(
        &self,
        target: &[u8],
        scratch: &mut RefindScratch,
        limit: usize,
    ) -> Option<(usize, usize, usize)> {
        let pattern = self
            .forward
            .get_or_init(|| GlobalPattern::new(self.query.iter().copied()));
        pattern.reset(&mut scratch.states);
        let states = &mut scratch.states;
        let blocks = states.len();
        let mut active = blocks;
        let mut bound = limit.min(self.length);
        let mut first = 0;
        let mut last = 0;
        for (index, &base) in target.iter().enumerate() {
            let offset = usize::from(pattern.slots[usize::from(base)]) * blocks;
            let equalities = &pattern.masks[offset..offset + blocks];
            // HW starts a new target substring for free at each column.
            let (mut positive_carry, mut negative_carry) = (0, 0);
            for ((state, &equal), &shift) in states[..active]
                .iter_mut()
                .zip(equalities)
                .zip(&pattern.shifts)
            {
                (positive_carry, negative_carry) =
                    advance_block_bits(state, equal, shift, positive_carry, negative_carry);
            }
            let carry = positive_carry as i8 - negative_carry as i8;
            if active < blocks
                && states[active - 1]
                    .distance
                    .wrapping_add_signed(-isize::from(carry))
                    <= bound
                && (equalities[active] & 1 != 0 || carry < 0)
            {
                // The next block can only re-enter through this boundary.
                let width = (self.length - active * 64).min(64);
                states[active] = GlobalState {
                    positive: u64::MAX,
                    negative: 0,
                    distance: (states[active - 1].distance + width)
                        .wrapping_add_signed(-isize::from(carry)),
                };
                advance_refind_block(
                    &mut states[active],
                    equalities[active],
                    pattern.boundaries[active],
                    carry,
                );
                active += 1;
            }
            if active == blocks {
                let distance = states[blocks - 1].distance;
                if distance <= bound {
                    if first == 0 || distance < bound {
                        first = index + 1;
                    }
                    last = index + 1;
                    bound = distance;
                }
            }
            // Only drop trailing blocks with no row at or below the bound.
            // Within a block, going up from its last row each row can only be
            // smaller than the one below by one where that row's vertical
            // delta is positive, so the block minimum is at least its last-row
            // distance minus the number of positive deltas (which is never
            // weaker than the 63-row bound `distance >= bound + 64`).
            // A last row at or below the bound keeps the block, so the
            // popcount is only needed when it is above.
            while active > 1 && {
                let state = &states[active - 1];
                let used = (pattern.boundaries[active - 1] << 1).wrapping_sub(1);
                state.distance > bound
                    && (state.distance >= bound + 64
                        || state
                            .distance
                            .saturating_sub((state.positive & used).count_ones() as usize)
                            > bound)
            } {
                active -= 1;
            }
        }
        (first != 0).then_some((bound, first, last))
    }

    fn recover(
        &self,
        target: &[u8],
        best: usize,
        first_end: usize,
        last_end: usize,
        scratch: &mut RefindScratch,
    ) -> RefindAlignment {
        if best == 0 {
            return RefindAlignment {
                distance: 0,
                start: last_end - self.length,
                end: last_end,
                terminal_insertions: 0,
            };
        }
        // Very narrow scalar bands cost less than visiting every 64-bit block.
        let blocks = (self.length - 1) / 64 + 1;
        if best.saturating_mul(2).saturating_add(1) < blocks.saturating_mul(4) {
            return self.recover_banded(target, best, first_end, last_end, scratch);
        }
        let forward = self
            .forward
            .get_or_init(|| GlobalPattern::new(self.query.iter().copied()));
        let reverse = self
            .reverse
            .get_or_init(|| GlobalPattern::new(self.query.iter().rev().copied()));
        let find_start = |end: usize, scratch: &mut RefindScratch| {
            reverse.reset(&mut scratch.states);
            let blocks = reverse.boundaries.len();
            // Only cells within `best` of the diagonal can be at or below
            // `best` (a cell is at least |i - j|), and only those decide
            // whether the last row equals `best`. Blocks entirely below that
            // band are started late from the upper-bound column `above + k`,
            // and blocks entirely above it are dropped with the top-row carry
            // of +1 per column fed to the next block. Both keep every cell an
            // upper bound on its true value, exact whenever the true value is
            // at most `best` (an optimal path never passes a larger cell), so
            // the columns where the last row equals `best` are unchanged.
            let mut low = 0;
            let mut high = blocks.min(best.saturating_add(1) / 64 + 1);
            let mut start = None;
            for j in 1..=end.min(self.length.saturating_add(best)) {
                while high < blocks && 64 * high < j + best {
                    let width = (self.length - high * 64).min(64);
                    scratch.states[high] = GlobalState {
                        positive: u64::MAX,
                        negative: 0,
                        distance: scratch.states[high - 1].distance + width,
                    };
                    high += 1;
                }
                while low + 1 < high && j > 64 * (low + 1) + best {
                    low += 1;
                }
                let distance =
                    reverse.advance_range(&mut scratch.states, target[end - j], low, high);
                if high == blocks && distance == best {
                    start = Some(end - j);
                }
            }
            start.expect("Optimal alignment has a start")
        };
        let start = find_start(last_end, scratch);
        let first_start = if first_end == last_end {
            start
        } else {
            find_start(first_end, scratch)
        };
        let first_target = &target[first_start..first_end];
        let width = first_target.len();
        forward.reset(&mut scratch.states);
        scratch.columns.clear();
        scratch.columns.reserve(width * forward.boundaries.len());
        // The traceback only visits cells at or below `best` and their
        // immediate neighbours, so the same band as `find_start` applies,
        // widened by a margin of two rows: blocks not yet started still hold
        // their column-zero values (an upper bound for rows below the band),
        // and dropped blocks are never read.
        let blocks = forward.boundaries.len();
        let mut low = 0;
        let mut high = blocks.min(best.saturating_add(3) / 64 + 1);
        for (column, &base) in first_target.iter().enumerate() {
            let j = column + 1;
            while high < blocks && 64 * high < j + best + 2 {
                let width = (self.length - high * 64).min(64);
                scratch.states[high] = GlobalState {
                    positive: u64::MAX,
                    negative: 0,
                    distance: scratch.states[high - 1].distance + width,
                };
                high += 1;
            }
            while low + 1 < high && j > 64 * (low + 1) + best + 2 {
                low += 1;
            }
            forward.advance_range(&mut scratch.states, base, low, high);
            scratch.columns.extend_from_slice(&scratch.states);
        }
        let (mut i, mut j) = (self.length, width);
        let mut trailing = 0;
        let mut leading = 0;
        let mut at_end = true;
        while i > 0 || j > 0 {
            let distance = scratch.distance(forward, i, j);
            if i > 0 && (j == 0 || distance == scratch.distance(forward, i - 1, j) + 1) {
                if at_end {
                    trailing += 1;
                }
                leading += 1;
                i -= 1;
            } else {
                at_end = false;
                leading = 0;
                if i == 0 || distance == scratch.distance(forward, i, j - 1) + 1 {
                    j -= 1;
                } else {
                    i -= 1;
                    j -= 1;
                }
            }
        }
        RefindAlignment {
            distance: best,
            start,
            end: last_end,
            terminal_insertions: trailing + leading,
        }
    }

    fn recover_banded(
        &self,
        target: &[u8],
        best: usize,
        first_end: usize,
        last_end: usize,
        scratch: &mut RefindScratch,
    ) -> RefindAlignment {
        let query = &self.query;
        let cap = best + 1;
        let find_start = |end: usize, scratch: &mut RefindScratch| {
            scratch.previous.clear();
            scratch
                .previous
                .extend((0..=self.length).map(|i| i.min(cap)));
            scratch.current.resize(self.length + 1, cap);
            scratch.current.fill(cap);
            let mut start = None;
            for j in 1..=end.min(self.length.saturating_add(best)) {
                scratch.current[0] = j.min(cap);
                let first = j.saturating_sub(best).max(1);
                let last = j.saturating_add(best).min(self.length);
                if first > 1 {
                    scratch.current[first - 1] = cap;
                }
                if last < self.length {
                    scratch.current[last + 1] = cap;
                }
                for i in first..=last {
                    let diagonal = scratch.previous[i - 1]
                        + usize::from(!refind_equal(query[self.length - i], target[end - j]));
                    scratch.current[i] = diagonal
                        .min(scratch.previous[i] + 1)
                        .min(scratch.current[i - 1] + 1)
                        .min(cap);
                }
                if scratch.current[self.length] == best {
                    start = Some(end - j);
                }
                std::mem::swap(&mut scratch.previous, &mut scratch.current);
            }
            start.expect("Optimal alignment has a start")
        };
        let start = find_start(last_end, scratch);
        let first_start = if first_end == last_end {
            start
        } else {
            find_start(first_end, scratch)
        };
        let first_target = &target[first_start..first_end];
        let width = first_target.len();
        scratch.previous.clear();
        scratch.previous.extend((0..=width).map(|j| j.min(cap)));
        scratch.current.resize(width + 1, cap);
        scratch.current.fill(cap);
        scratch.paths.clear();
        scratch.row_offsets.resize(self.length + 1, 0);
        scratch.row_starts.resize(self.length + 1, 0);
        for i in 1..=self.length {
            scratch.current[0] = i.min(cap);
            let first = i.saturating_sub(best).max(1);
            let last = i.saturating_add(best).min(width);
            if first > 1 && first - 1 <= width {
                scratch.current[first - 1] = cap;
            }
            if last < width {
                scratch.current[last + 1] = cap;
            }
            scratch.row_offsets[i] = scratch.paths.len();
            scratch.row_starts[i] = first;
            for j in first..=last {
                let insertion = scratch.previous[j] + 1;
                let deletion = scratch.current[j - 1] + 1;
                let diagonal = scratch.previous[j - 1]
                    + usize::from(!refind_equal(query[i - 1], first_target[j - 1]));
                scratch.current[j] = insertion.min(deletion).min(diagonal).min(cap);
                scratch.paths.push(if scratch.current[j] == insertion {
                    1
                } else if scratch.current[j] == deletion {
                    2
                } else {
                    0
                });
            }
            std::mem::swap(&mut scratch.previous, &mut scratch.current);
        }
        let (mut i, mut j) = (self.length, width);
        let mut trailing = 0;
        let mut leading = 0;
        let mut at_end = true;
        while i > 0 || j > 0 {
            let op = if i == 0 {
                2
            } else if j == 0 {
                1
            } else {
                scratch.paths[scratch.row_offsets[i] + j - scratch.row_starts[i]]
            };
            if op == 1 {
                if at_end {
                    trailing += 1;
                }
                leading += 1;
                i -= 1;
            } else {
                at_end = false;
                leading = 0;
                if op == 0 {
                    i -= 1;
                }
                j -= 1;
            }
        }
        RefindAlignment {
            distance: best,
            start,
            end: last_end,
            terminal_insertions: trailing + leading,
        }
    }
}

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

    fn random_sequence(length: usize, alphabet: &[u8], state: &mut u64) -> Vec<u8> {
        (0..length)
            .map(|_| {
                *state ^= *state << 13;
                *state ^= *state >> 7;
                *state ^= *state << 17;
                alphabet[*state as usize % alphabet.len()]
            })
            .collect()
    }

    #[test]
    fn limited_single_pass_is_exact_or_reports_nothing() {
        // The limited scan must return exactly the unlimited endpoints whenever
        // the optimum is within the limit, and nothing otherwise.
        let mut seed = 0x51ed_2701_9c3b_77a1u64;
        let mut scratch = RefindScratch::default();
        for case in 0..900 {
            let length = [
                1, 5, 40, 63, 64, 65, 100, 127, 128, 129, 200, 300, 520, 1100,
            ][case % 14];
            let alphabet: &[u8] = if case % 3 == 0 { b"ACGTNEX" } else { b"ACGT" };
            let query = random_sequence(length, alphabet, &mut seed);
            let mut target = random_sequence(length * 4 + 31, alphabet, &mut seed);
            // Plant a mutated copy of the query (substitutions and indels).
            if case % 5 != 0 {
                let mut copy = query.clone();
                let rate = [0usize, 50, 20, 9, 4][case % 5];
                if rate > 0 {
                    let mut index = (seed as usize) % rate.max(1);
                    while index < copy.len() {
                        seed ^= seed << 13;
                        seed ^= seed >> 7;
                        seed ^= seed << 17;
                        match seed % 3 {
                            0 => copy[index] = b"ACGT"[(seed >> 8) as usize % 4],
                            1 => {
                                copy.insert(index, b"ACGT"[(seed >> 8) as usize % 4]);
                            }
                            _ => {
                                copy.remove(index);
                            }
                        }
                        index += 1 + (seed >> 16) as usize % rate;
                    }
                }
                let at = (seed as usize) % (target.len() - copy.len().min(target.len()) + 1);
                let end = (at + copy.len()).min(target.len());
                target.splice(at..end, copy.iter().copied());
                if case % 7 == 0 {
                    // A second, identical copy creates endpoint ties.
                    target.extend(&copy);
                }
            }
            if case % 4 == 0 {
                let padding = length / 2;
                target = [vec![b'E'; padding], target, vec![b'E'; padding]].concat();
            }
            let matcher = RefindMatcher::new(&query).unwrap();
            let full = matcher
                .single_pass_endpoints_within(&target, &mut scratch, length)
                .unwrap();
            let mut expected = (usize::MAX, 0, 0);
            for (end, distance) in matcher.matcher.find_all_end(&target, length) {
                if distance < expected.0 {
                    expected = (distance, end + 1, end + 1);
                } else if distance == expected.0 {
                    expected.2 = end + 1;
                }
            }
            assert_eq!(full, expected, "case {case}");
            assert_eq!(
                matcher.single_pass_endpoints(&target, &mut scratch),
                Some(full)
            );
            let limits = [0, 1, 2, 3, 8, 16, 33, 63, 64, 65, 100, 130, 200, 400];
            for limit in limits.into_iter().chain([
                full.0.saturating_sub(1),
                full.0,
                full.0 + 1,
                length / 8 + 8,
            ]) {
                let limited = matcher.single_pass_endpoints_within(&target, &mut scratch, limit);
                if full.0 <= limit {
                    assert_eq!(limited, Some(full), "case {case} limit {limit}");
                } else {
                    assert_eq!(limited, None, "case {case} limit {limit}");
                }
            }
        }
    }

    #[test]
    fn single_pass_matches_full_scan_across_blocks_and_ties() {
        let mut seed = 0x2d98_039f_b369_a3acu64;
        let mut scratch = RefindScratch::default();
        for case in 0..1800 {
            let length = [1, 2, 63, 64, 65, 127, 128, 129, 256, 511, 1024][case % 11];
            let alphabet: &[u8] = if case % 2 == 0 { b"ACGT" } else { b"ACGTNEX" };
            let query = random_sequence(length, alphabet, &mut seed);
            let mut target = random_sequence(length * 3 + 7, alphabet, &mut seed);
            match case % 7 {
                0 => target.truncate(1),
                1 => {
                    target[..length].copy_from_slice(&query);
                    target[2 * length..3 * length].copy_from_slice(&query);
                }
                2 => {
                    target[..length].copy_from_slice(&query);
                    target[2 * length..3 * length].copy_from_slice(&query);
                    target[length / 2] = b'Z';
                    target[2 * length + length / 2] = b'Z';
                }
                3 => drop(target.splice(0..0, vec![b'E'; length / 2])),
                4 => target.extend(&query),
                _ => {}
            }
            let matcher = RefindMatcher::new(&query).unwrap();
            let mut expected = (usize::MAX, 0, 0);
            for (end, distance) in matcher.matcher.find_all_end(&target, length) {
                if distance < expected.0 {
                    expected = (distance, end + 1, end + 1);
                } else if distance == expected.0 {
                    expected.2 = end + 1;
                }
            }
            assert_eq!(
                matcher.single_pass_endpoints(&target, &mut scratch),
                Some(expected),
                "case {case}"
            );
            if case < 220 {
                assert_eq!(
                    matcher.align_with_scratch(&target, &mut scratch),
                    Some(scalar_recovery_comparator(&matcher, &target)),
                    "recovered alignment for case {case}"
                );
            }
        }
    }

    #[test]
    fn compact_global_masks_preserve_matches_for_every_target_byte() {
        let all_bytes: Vec<_> = (u8::MIN..=u8::MAX).collect();
        let mut sparse_query = b"ACGTNEX".repeat(10);
        sparse_query.push(b'A');
        for query in [&sparse_query, &all_bytes] {
            for query in [query.clone(), query.iter().rev().copied().collect()] {
                let pattern = GlobalPattern::new(query.iter().copied());
                let blocks = query.len().div_ceil(64);
                for target in u8::MIN..=u8::MAX {
                    let slot = usize::from(pattern.slots[usize::from(target)]);
                    for (block, bases) in query.chunks(64).enumerate() {
                        let expected = bases.iter().enumerate().fold(0, |mask, (i, &base)| {
                            mask | (u64::from(refind_equal(base, target)) << i)
                        });
                        assert_eq!(pattern.masks[slot * blocks + block], expected);
                    }
                }
                if query.len() == sparse_query.len() {
                    assert_eq!(pattern.masks.len(), 8 * blocks);
                }
            }
        }
    }

    #[test]
    fn global_bit_columns_match_scalar_cells_across_block_boundaries() {
        let mut seed = 0x73ac_49b1_02u64;
        let mut scratch = RefindScratch::default();
        for length in [1, 2, 63, 64, 65, 127, 128, 129, 257] {
            let query = random_sequence(length, b"ACGTNEXac?", &mut seed);
            let target = random_sequence(length + 17, b"ACGTNEXac?", &mut seed);
            let pattern = GlobalPattern::new(query.iter().copied());
            pattern.reset(&mut scratch.states);
            scratch.columns.clear();
            let mut previous: Vec<_> = (0..=length).collect();
            let mut current = vec![0; length + 1];
            for (column, &base) in target.iter().enumerate() {
                pattern.advance(&mut scratch.states, base);
                scratch.columns.extend_from_slice(&scratch.states);
                current[0] = column + 1;
                for row in 1..=length {
                    current[row] = (previous[row] + 1)
                        .min(current[row - 1] + 1)
                        .min(previous[row - 1] + usize::from(!refind_equal(query[row - 1], base)));
                }
                for (row, &expected) in current.iter().enumerate() {
                    assert_eq!(
                        scratch.distance(&pattern, row, column + 1),
                        expected,
                        "length={length} row={row} column={column}"
                    );
                }
                std::mem::swap(&mut previous, &mut current);
            }
        }
    }

    #[test]
    fn adaptive_refind_preserves_ambiguous_indels_ties_and_reused_state() {
        let mut seed = 0x8bc4_5921_dfe3u64;
        let mut scratch = RefindScratch::default();
        for case in 0..400 {
            let length = 1 + case % 260;
            let query = random_sequence(length, b"ACGTNEX", &mut seed);
            let mut target = match case % 4 {
                0 => random_sequence(length * 2 + 9, b"ACGTNEX", &mut seed),
                1 => query.clone(),
                2 => query[..length / 2 + 1].to_vec(),
                _ => [query.as_slice(), b"TTGAC", query.as_slice()].concat(),
            };
            if case % 4 != 0 {
                for index in (0..target.len()).step_by(7) {
                    target[index] = b"ACGT"[index % 4];
                }
                target.insert(target.len() / 2, b'C');
                target.remove(target.len() / 3);
            }
            let padding = length / 2;
            let target = [vec![b'E'; padding], target, vec![b'E'; padding]].concat();
            let matcher = RefindMatcher::new(&query).unwrap();
            let expected = full_refind(&query, &target);
            assert_eq!(
                matcher.align_with_scratch(&target, &mut scratch),
                expected,
                "case={case} query={query:?} target={target:?}"
            );
            assert_eq!(matcher.align_with_scratch(b"", &mut scratch), None);
        }
        for length in [63, 64, 65, 127, 128, 129, 1024] {
            let query = random_sequence(length, b"ACGT", &mut seed);
            let mut target = random_sequence(length * 3, b"ACGT", &mut seed);
            target.extend(&query);
            target.extend(&query);
            let matcher = RefindMatcher::new(&query).unwrap();
            assert_eq!(
                matcher.align_with_scratch(&target, &mut scratch),
                full_refind(&query, &target)
            );
        }
    }

    #[test]
    fn block_recovery_band_matches_full_traceback() {
        // Multi-block queries with divergence high enough that recovery takes
        // the bit-parallel (not the scalar banded) path, so the band limits in
        // `find_start` and the forward pass are exercised at many offsets,
        // including copies near both ends of the padded target and endpoint ties.
        let mut seed = 0x6a09_e667_f3bc_c908u64;
        let mut scratch = RefindScratch::default();
        for case in 0..360 {
            let length = [65, 100, 127, 128, 129, 191, 192, 193, 300, 520, 700][case % 11];
            let alphabet: &[u8] = if case % 3 == 0 { b"ACGTNEX" } else { b"ACGT" };
            let query = random_sequence(length, alphabet, &mut seed);
            let mut copy = query.clone();
            let rate = [2usize, 3, 4, 6, 9, 15][case % 6];
            let mut index = 0;
            while index < copy.len() {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                match seed % 3 {
                    0 => copy[index] = b"ACGT"[(seed >> 8) as usize % 4],
                    1 => copy.insert(index, b"ACGT"[(seed >> 8) as usize % 4]),
                    _ => drop(copy.remove(index)),
                }
                index += 1 + (seed >> 16) as usize % rate;
            }
            let flank = [0, 3, length / 3, length][case % 4];
            let mut target = random_sequence(flank, b"ACGT", &mut seed);
            target.extend(&copy);
            target.extend(random_sequence([length, 0, 7][case % 3], b"ACGT", &mut seed));
            if case % 5 == 0 {
                target.extend(&copy);
            }
            let padding = length / 2;
            let target = [vec![b'E'; padding], target, vec![b'E'; padding]].concat();
            let matcher = RefindMatcher::new(&query).unwrap();
            assert_eq!(
                matcher.align_with_scratch(&target, &mut scratch),
                full_refind(&query, &target),
                "case={case} length={length} rate={rate}"
            );
        }
    }

    fn scalar_recovery_comparator(matcher: &RefindMatcher, target: &[u8]) -> RefindAlignment {
        let mut best = usize::MAX;
        let mut first_end = 0;
        let mut last_end = 0;
        for (end, distance) in matcher.matcher.find_all_end(target, matcher.length) {
            if distance < best {
                best = distance;
                first_end = end + 1;
                last_end = end + 1;
            } else if distance == best {
                last_end = end + 1;
            }
        }
        if best == 0 {
            return RefindAlignment {
                distance: 0,
                start: last_end - matcher.length,
                end: last_end,
                terminal_insertions: 0,
            };
        }
        matcher.recover_banded(
            target,
            best,
            first_end,
            last_end,
            &mut RefindScratch::default(),
        )
    }

    #[test]
    #[ignore = "manual release-mode synthetic timing; not a whole-pipeline benchmark"]
    fn refind_kernel_timing() {
        use std::hint::black_box;
        use std::time::{Duration, Instant};
        let mut seed = 0x783c_1bfe_943du64;
        for (label, length, window, change_every) in [
            ("near-match", 1200, 12000, 50),
            ("divergent-match", 1200, 12000, 5),
            ("unrelated", 1200, 12000, 0),
            ("long-divergent", 3000, 16000, 5),
        ] {
            let query = random_sequence(length, b"ACGT", &mut seed);
            let mut target = random_sequence(window, b"ACGT", &mut seed);
            if change_every != 0 {
                let offset = window / 3;
                target[offset..offset + length].copy_from_slice(&query);
                for position in (offset..offset + length).step_by(change_every) {
                    target[position] = match target[position] {
                        b'A' => b'C',
                        b'C' => b'G',
                        b'G' => b'T',
                        _ => b'A',
                    };
                }
            }
            let target = [vec![b'E'; length / 2], target, vec![b'E'; length / 2]].concat();
            let matcher = RefindMatcher::new(&query).unwrap();
            let mut scratch = RefindScratch::default();
            assert_eq!(
                matcher.align_with_scratch(&target, &mut scratch).unwrap(),
                scalar_recovery_comparator(&matcher, &target)
            );
            let mut scalar = Vec::new();
            let mut optimized = Vec::new();
            for round in 0..7 {
                for old in if round % 2 == 0 {
                    [true, false]
                } else {
                    [false, true]
                } {
                    let start = Instant::now();
                    for _ in 0..4 {
                        if old {
                            black_box(scalar_recovery_comparator(
                                black_box(&matcher),
                                black_box(&target),
                            ));
                        } else {
                            black_box(matcher.align_with_scratch(black_box(&target), &mut scratch));
                        }
                    }
                    if old {
                        scalar.push(start.elapsed());
                    } else {
                        optimized.push(start.elapsed());
                    }
                }
            }
            scalar.sort();
            optimized.sort();
            let per_call = |duration: Duration| duration.as_secs_f64() * 1000.0 / 4.0;
            eprintln!(
                "{label}: scalar_recovery_ms={:.3} optimized_ms={:.3} ratio={:.2}x",
                per_call(scalar[3]),
                per_call(optimized[3]),
                scalar[3].as_secs_f64() / optimized[3].as_secs_f64()
            );
        }
    }

    fn full_refind(query: &[u8], target: &[u8]) -> Option<RefindAlignment> {
        if query.is_empty() || target.is_empty() {
            return None;
        }
        let matcher = RefindMatcher::new(query)?;
        let mut best = usize::MAX;
        let mut first_end = 0;
        let mut last_end = 0;
        for (end, distance) in matcher.matcher.find_all_end(target, query.len()) {
            if distance < best {
                best = distance;
                first_end = end + 1;
                last_end = end + 1;
            } else if distance == best {
                last_end = end + 1;
            }
        }
        if best == usize::MAX {
            return None;
        }
        if best == 0 {
            return Some(RefindAlignment {
                distance: 0,
                start: last_end - query.len(),
                end: last_end,
                terminal_insertions: 0,
            });
        }
        let find_start = |end: usize| {
            let cap = best + 1;
            let mut previous: Vec<usize> = (0..=query.len()).map(|i| i.min(cap)).collect();
            let mut current = vec![cap; query.len() + 1];
            let mut start = None;
            for j in 1..=end.min(query.len().saturating_add(best)) {
                current.fill(cap);
                current[0] = j.min(cap);
                for i in j.saturating_sub(best).max(1)..=j.saturating_add(best).min(query.len()) {
                    let diagonal = previous[i - 1]
                        + usize::from(!refind_equal(query[query.len() - i], target[end - j]));
                    current[i] = diagonal
                        .min(previous[i] + 1)
                        .min(current[i - 1] + 1)
                        .min(cap);
                }
                if current[query.len()] == best {
                    start = Some(end - j);
                }
                std::mem::swap(&mut previous, &mut current);
            }
            start.expect("optimal alignment has a start")
        };
        let start = find_start(last_end);
        let first_start = if first_end == last_end {
            start
        } else {
            find_start(first_end)
        };
        let first_target = &target[first_start..first_end];
        let width = first_target.len();
        let cap = best + 1;
        let mut previous: Vec<usize> = (0..=width).map(|j| j.min(cap)).collect();
        let mut current = vec![cap; width + 1];
        let mut paths = vec![0u8; (query.len() + 1) * (width + 1)];
        for i in 1..=query.len() {
            current.fill(cap);
            current[0] = i.min(cap);
            paths[i * (width + 1)] = 1;
            for j in i.saturating_sub(best).max(1)..=i.saturating_add(best).min(width) {
                let insertion = previous[j] + 1;
                let deletion = current[j - 1] + 1;
                let diagonal =
                    previous[j - 1] + usize::from(!refind_equal(query[i - 1], first_target[j - 1]));
                current[j] = insertion.min(deletion).min(diagonal).min(cap);
                paths[i * (width + 1) + j] = if current[j] == insertion {
                    1
                } else if current[j] == deletion {
                    2
                } else {
                    0
                };
            }
            std::mem::swap(&mut previous, &mut current);
        }
        let (mut i, mut j) = (query.len(), width);
        let mut operations = Vec::new();
        while i > 0 || j > 0 {
            let op = if i == 0 {
                2
            } else {
                paths[i * (width + 1) + j]
            };
            operations.push(op);
            match op {
                0 => {
                    i -= 1;
                    j -= 1;
                }
                1 => i -= 1,
                _ => j -= 1,
            }
        }
        Some(RefindAlignment {
            distance: best,
            start,
            end: last_end,
            terminal_insertions: operations.iter().take_while(|&&op| op == 1).count()
                + operations.iter().rev().take_while(|&&op| op == 1).count(),
        })
    }

    #[test]
    fn banded_refind_matches_full_traceback() {
        let mut sequences = vec![Vec::new()];
        for length in 1..=4 {
            for mut code in 0..3_usize.pow(length) {
                let mut sequence = Vec::new();
                for _ in 0..length {
                    sequence.push(b"ACN"[code % 3]);
                    code /= 3;
                }
                sequences.push(sequence);
            }
        }
        for query in &sequences {
            let matcher = RefindMatcher::new(query);
            for target in &sequences {
                let actual = matcher.as_ref().and_then(|matcher| matcher.align(target));
                assert_eq!(
                    actual,
                    full_refind(query, target),
                    "query={query:?}, target={target:?}"
                );
            }
        }
        for length in [63, 64, 65, 127, 128, 129, 257] {
            let query: Vec<_> = (0..length)
                .map(|i| b"ACGTNE"[(i * 7 + i / 3) % 6])
                .collect();
            let mut target = b"GGTTA".to_vec();
            target.extend(&query);
            target[length / 2] = b'C';
            target.remove(length / 3);
            target.insert(length / 4, b'T');
            target.extend(b"CCA");
            let matcher = RefindMatcher::new(&query).unwrap();
            assert_eq!(matcher.align(&target), full_refind(&query, &target));
        }
    }

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
