//! Read-only reference evidence helpers ported from Rustaroo.

use super::align::{global_identity, terminal_identity_parts};
use super::fxhash::{HashMap, HashSet};
use super::genetic_code::orf_status;
use super::model::{Call, Pools};
use super::reference::{join_fragments, Allele};

const ALLELE_EDGE_BP: usize = 90;
const ABUT_GAP: i32 = 3;
const TERMINAL_IDENTITY: f64 = 0.98;
pub(super) fn alleles_of(
    indexes: &[usize],
    calls: &[Call],
    pools: &Pools,
    contig_lengths: &[Vec<usize>],
    locus_ids: &[u32],
) -> Vec<Allele> {
    // Two abutting pieces of one locus are a broken gene, not a shorter
    // finished stop. Counting either piece would let the break define the
    // reference and throw away the intact copy.
    let mut parts_of_locus: HashMap<u32, u32> = HashMap::default();
    for &index in indexes {
        *parts_of_locus.entry(locus_ids[index]).or_default() += 1;
    }
    let mut samples_for_dna: HashMap<u32, HashSet<u32>> = HashMap::default();
    let mut finished_for_dna: HashMap<u32, HashSet<u32>> = HashMap::default();
    let mut closed: HashMap<u32, bool> = HashMap::default();
    for &index in indexes {
        let call = &calls[index];
        samples_for_dna
            .entry(call.dna)
            .or_default()
            .insert(call.sample);
        let closed = *closed.entry(call.dna).or_insert_with(|| {
            let status = orf_status(&pools.dna[call.dna as usize]);
            status.is_intact() && status.terminal_stop
        });
        let whole = parts_of_locus.get(&locus_ids[index]).copied().unwrap_or(0) == 1;
        // One internal whole copy is enough for that sample to vouch for the
        // end. A second copy parked on a contig break does not cancel it.
        if closed && whole && !call.refound && !contig_clipped(call, contig_lengths) {
            finished_for_dna
                .entry(call.dna)
                .or_default()
                .insert(call.sample);
        }
    }
    let mut alleles: Vec<Allele> = samples_for_dna
        .into_iter()
        .map(|(dna, samples)| Allele {
            dna: pools.dna[dna as usize].clone(),
            count: samples.len() as u32,
            synthetic: false,
            finished: finished_for_dna
                .get(&dna)
                .map(|samples| samples.len() as u32)
                .unwrap_or(0),
        })
        .collect();
    alleles.sort_by(|a, b| b.count.cmp(&a.count).then(a.dna.cmp(&b.dna)));
    alleles
}

/// A copy this close to either contig boundary may have been clipped. It
/// stays in the family, and it does not get to define where the gene ends.
fn contig_clipped(call: &Call, contig_lengths: &[Vec<usize>]) -> bool {
    let Some(length) = contig_lengths
        .get(call.sample as usize)
        .and_then(|contigs| contigs.get(call.contig as usize))
        .copied()
    else {
        return true;
    };
    // FIX(orientation): both sides count the bases outside the call
    // (`start - 1` before it, `length - end` after it). The old `start <= 90`
    // test was one base stricter on the left than on the right, so a copy
    // exactly 90 bp from an end was clipped in one orientation only.
    (call.start as usize).saturating_sub(1) <= ALLELE_EDGE_BP
        || length.saturating_sub(call.end as usize) <= ALLELE_EDGE_BP
}

pub(super) fn joined_allele(
    indexes: &[usize],
    calls: &[Call],
    pools: &Pools,
    locus_ids: &[u32],
) -> Option<String> {
    let mut joins: HashMap<String, u32> = HashMap::default();
    for (joined, _) in joined_models(indexes, calls, pools, locus_ids) {
        *joins.entry(joined).or_default() += 1;
    }
    joins
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
        .map(|(dna, _)| dna)
}

pub(super) fn joined_models(
    indexes: &[usize],
    calls: &[Call],
    pools: &Pools,
    locus_ids: &[u32],
) -> Vec<(String, Vec<usize>)> {
    let mut by_locus: HashMap<(u32, u32, u32), Vec<usize>> = HashMap::default();
    for &index in indexes {
        let call = &calls[index];
        by_locus
            .entry((call.sample, call.contig, locus_ids[index]))
            .or_default()
            .push(index);
    }
    let mut joins = Vec::new();
    for group in by_locus.values() {
        if group.len() < 2 {
            continue;
        }
        let mut ordered = group.clone();
        ordered.sort_by_key(|index| calls[*index].start);
        let strand = calls[ordered[0]].strand;
        if ordered.iter().any(|index| calls[*index].strand != strand) {
            continue;
        }
        if strand < 0 {
            ordered.reverse();
        }
        let mut abutting = true;
        for pair in ordered.windows(2) {
            let gap = if strand >= 0 {
                i64::from(calls[pair[1]].start) - i64::from(calls[pair[0]].end) - 1
            } else {
                i64::from(calls[pair[0]].start) - i64::from(calls[pair[1]].end) - 1
            };
            if !(0..=3).contains(&gap) {
                abutting = false;
                break;
            }
        }
        if !abutting {
            continue;
        }
        let parts: Vec<&str> = ordered
            .iter()
            .map(|index| pools.dna[calls[*index].dna as usize].as_str())
            .collect();
        if let Some(joined) = join_fragments(&parts) {
            joins.push((joined, ordered));
        }
    }
    joins
}

pub(super) fn cluster_reps(
    groups: &HashMap<u32, Vec<usize>>,
    calls: &[Call],
    pools: &Pools,
) -> HashMap<u32, Vec<u8>> {
    let mut reps = HashMap::default();
    for (&cluster, indexes) in groups {
        let protein = indexes
            .iter()
            .map(|index| &pools.protein[calls[*index].protein as usize])
            // FIX(orientation): longest, then smallest sequence on a length tie.
            .max_by(|a, b| a.len().cmp(&b.len()).then_with(|| b.cmp(a)))
            .cloned()
            .unwrap_or_default();
        reps.insert(cluster, protein);
    }
    reps
}

pub(super) fn membership_loci(calls: &[Call], pools: &Pools) -> Vec<u32> {
    let mut ids = vec![0; calls.len()];
    for (id, locus) in physical_loci(calls, pools).into_iter().enumerate() {
        ids[locus.first..=locus.last].fill(id as u32);
    }
    ids
}

struct Locus {
    first: usize,
    last: usize,
    end: u32,
}

fn same_physical_locus(
    left: usize,
    right: usize,
    calls: &[Call],
    pools: &Pools,
    template: &[u8],
    cache: &mut HashMap<(u32, u32), bool>,
) -> bool {
    let a = &calls[left];
    let b = &calls[right];
    if right != left + 1 || a.sample != b.sample || a.contig != b.contig || a.strand != b.strand {
        return false;
    }
    if a.locus.is_some() && a.locus == b.locus {
        return true;
    }
    let overlap = i64::from(a.end.min(b.end)) - i64::from(a.start.max(b.start)) + 1;
    let shorter =
        (i64::from(a.end) - i64::from(a.start) + 1).min(i64::from(b.end) - i64::from(b.start) + 1);
    if shorter > 0 && overlap * 100 >= shorter * 95 {
        return true;
    }
    let gap = i64::from(b.start) - i64::from(a.end) - 1;
    if gap.abs() > i64::from(ABUT_GAP) {
        return false;
    }
    let a_protein = &pools.protein[a.protein as usize];
    let b_protein = &pools.protein[b.protein as usize];
    let (head, tail) = if a.strand >= 0 {
        (a_protein.as_slice(), b_protein.as_slice())
    } else {
        (b_protein.as_slice(), a_protein.as_slice())
    };
    // Repetitive fragments can resemble each other. A longer observed allele
    // can prove their join, whereas two full tandem copies exceed its length.
    if template.len() > head.len().max(tail.len())
        && (head.len() + tail.len()) as f64 >= 0.95 * template.len() as f64
        && terminal_identity_parts(head, tail, template)
            .is_some_and(|identity| identity + f64::EPSILON >= TERMINAL_IDENTITY)
    {
        return true;
    }
    if a.protein == b.protein {
        return false;
    }
    let pair = (a.protein.min(b.protein), a.protein.max(b.protein));
    *cache.entry(pair).or_insert_with(|| {
        global_identity(
            &pools.protein[a.protein as usize],
            &pools.protein[b.protein as usize],
        )
        .is_some_and(|identity| identity < 0.4)
    })
}

fn physical_loci(calls: &[Call], pools: &Pools) -> Vec<Locus> {
    let representatives = representative_ids(calls, pools);
    let mut loci: Vec<Locus> = Vec::new();
    let mut cache = HashMap::default();
    for (index, call) in calls.iter().enumerate() {
        if let Some(previous) = loci.last_mut() {
            if calls[previous.first].cluster == call.cluster
                && same_physical_locus(
                    previous.last,
                    index,
                    calls,
                    pools,
                    &pools.protein[representatives[&call.cluster] as usize],
                    &mut cache,
                )
            {
                previous.last = index;
                previous.end = previous.end.max(call.end);
                continue;
            }
        }
        loci.push(Locus {
            first: index,
            last: index,
            end: call.end,
        });
    }
    loci
}

fn representative_ids(calls: &[Call], pools: &Pools) -> HashMap<u32, u32> {
    let mut best = HashMap::default();
    for call in calls {
        let slot = best.entry(call.cluster).or_insert(call.protein);
        // FIX(orientation): equal lengths go to the smaller sequence, not to
        // whichever copy came first in coordinate order.
        let (candidate, current) = (
            &pools.protein[call.protein as usize],
            &pools.protein[*slot as usize],
        );
        if candidate.len() > current.len()
            || (candidate.len() == current.len() && candidate < current)
        {
            *slot = call.protein;
        }
    }
    best
}
