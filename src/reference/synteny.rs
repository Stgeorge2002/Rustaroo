//! Sequence-only tail screening from Rustaroo; no family or graph construction.

use super::fxhash::HashMap;

pub(crate) struct FusionIndex {
    postings: HashMap<[u8; 4], Vec<u32>>,
    tiny: Vec<u32>,
}

impl FusionIndex {
    pub(crate) fn new(reps: &HashMap<u32, Vec<u8>>, cluster_ids: &[u32]) -> Self {
        let mut postings: HashMap<[u8; 4], Vec<u32>> = HashMap::default();
        let mut tiny = Vec::new();
        for &id in cluster_ids {
            let protein = &reps[&id];
            if protein.len() < 4 {
                tiny.push(id);
                continue;
            }
            for w in protein.windows(4) {
                let list = postings.entry([w[0], w[1], w[2], w[3]]).or_default();
                if list.last() != Some(&id) {
                    list.push(id);
                }
            }
        }
        Self { postings, tiny }
    }

    pub(crate) fn sharing(&self, sequence: &[u8], include_tiny: bool) -> Vec<u32> {
        let mut out = Vec::new();
        if sequence.len() >= 4 {
            for w in sequence.windows(4) {
                if let Some(list) = self.postings.get(&[w[0], w[1], w[2], w[3]]) {
                    out.extend_from_slice(list);
                }
            }
        }
        if include_tiny {
            out.extend_from_slice(&self.tiny);
        }
        out.sort_unstable();
        out.dedup();
        out
    }
}

pub(crate) type Composition = [u32; 32];

pub(crate) fn composition(sequence: &[u8]) -> Composition {
    let mut counts = [0u32; 32];
    for &residue in sequence {
        counts[(residue & 31) as usize] += 1;
    }
    counts
}

pub(crate) fn composition_bound(a: &Composition, b: &Composition) -> usize {
    a.iter().zip(b).map(|(x, y)| (*x).min(*y) as usize).sum()
}
