//! Iterative CD-HIT clustering at multiple thresholds
//!

use crate::alignment::HwEditScratch;
use crate::cdhit;
use crate::network::Graph;
use crate::network::{GeneEdge, GeneNode};
use crate::single_linkage::single_linkage_clustering_by;
use anyhow::Result;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;
use tempfile::TempDir;

/// Run iterative CD-HIT at multiple thresholds
pub fn iterative_cdhit(
    graph: &Graph<GeneNode, GeneEdge>,
    nodes: &[NodeIndex],
    outdir: &str,
    thresholds: &[f64],
    len_dif_percent: f64,
    n_cpu: usize,
    _quiet: bool,
    dna: bool,
) -> Result<Vec<Vec<usize>>> {
    if nodes.is_empty() {
        return Ok(vec![]);
    }

    let temp_dir = TempDir::new_in(outdir)?;
    let temp_input = temp_dir.path().join("input.fasta");

    // Write sequences to temp file
    let mut sequence_count = 0;
    {
        let mut file = File::create(&temp_input)?;
        for (i, &node_idx) in nodes.iter().enumerate() {
            let node = &graph[node_idx];
            let seq = if dna {
                node.dna.get(node.max_len_id)
            } else {
                None
            }
            .or_else(|| node.protein.get(node.max_len_id));

            let Some(seq) = seq.filter(|sequence| !sequence.is_empty()) else {
                continue;
            };
            writeln!(file, ">{}", i)?;
            writeln!(file, "{}", seq)?;
            sequence_count += 1;
        }
    }

    if sequence_count == 0 {
        anyhow::bail!(
            "Cannot run iterative CD-HIT: none of the {} graph nodes have a usable sequence",
            nodes.len()
        );
    }

    // Run CD-HIT with each threshold
    let mut current_input = temp_input.to_str().unwrap().to_string();

    for (i, &threshold) in thresholds.iter().enumerate() {
        let output = temp_dir.path().join(format!("output_{}.fasta", i));
        let output_str = output.to_str().unwrap();

        if dna {
            cdhit::run_cdhit_est(
                &current_input,
                output_str,
                threshold,
                len_dif_percent,
                1, // strand: both
                n_cpu,
            )?;
        } else {
            cdhit::run_cdhit(
                &current_input,
                output_str,
                threshold,
                len_dif_percent,
                n_cpu,
                outdir,
            )?;
        }

        current_input = output_str.to_string();
    }

    // Parse final clusters
    let clstr_file = format!("{}.clstr", current_input);
    let clusters = cdhit::parse_cdhit_clusters(&clstr_file)?;

    // Convert to node index clusters
    let mut result_clusters = Vec::new();
    for (_cluster_id, members) in clusters {
        let node_indices: Vec<usize> = members
            .iter()
            .filter_map(|member_str| {
                member_str.parse::<usize>().ok().and_then(|idx| {
                    if idx < nodes.len() {
                        Some(nodes[idx].index())
                    } else {
                        None
                    }
                })
            })
            .collect();

        if !node_indices.is_empty() {
            result_clusters.push(node_indices);
        }
    }

    Ok(result_clusters)
}

const MAX_CACHED_SEQUENCES: usize = 65_536;
const MAX_CACHED_PAIRS: usize = 250_000;
const QGRAM_LENGTH: usize = 3;

#[derive(Default)]
pub struct PairwiseDistanceCache {
    sequence_ids: HashMap<Arc<str>, u32>,
    sequence_qgrams: Vec<Vec<u32>>,
    distances: HashMap<(u32, u32, usize), Option<usize>>,
    scratch: HwEditScratch,
    hits: usize,
    calculations: usize,
    qgram_rejections: usize,
}

impl PairwiseDistanceCache {
    pub fn hit_count(&self) -> usize {
        self.hits
    }

    pub fn calculation_count(&self) -> usize {
        self.calculations
    }

    pub fn qgram_rejection_count(&self) -> usize {
        self.qgram_rejections
    }

    fn sequence_id(&mut self, sequence: &str) -> Option<u32> {
        if let Some(&id) = self.sequence_ids.get(sequence) {
            return Some(id);
        }
        if self.sequence_ids.len() >= MAX_CACHED_SEQUENCES {
            return None;
        }

        let id = self.sequence_ids.len() as u32;
        self.sequence_ids.insert(Arc::from(sequence), id);
        self.sequence_qgrams
            .push(encode_qgrams(sequence.as_bytes()));
        Some(id)
    }

    fn distance(
        &mut self,
        first_id: Option<u32>,
        second_id: Option<u32>,
        first: &str,
        second: &str,
        max_distance: usize,
    ) -> Option<usize> {
        if let (Some(first_id), Some(second_id)) = (first_id, second_id) {
            let key = (first_id, second_id, max_distance);
            if let Some(&distance) = self.distances.get(&key) {
                self.hits += 1;
                return distance;
            }

            if qgram_excess(
                &self.sequence_qgrams[first_id as usize],
                &self.sequence_qgrams[second_id as usize],
            ) > QGRAM_LENGTH.saturating_mul(max_distance)
            {
                self.qgram_rejections += 1;
                if self.distances.len() < MAX_CACHED_PAIRS {
                    self.distances.insert(key, None);
                }
                return None;
            }
            if first == second {
                if self.distances.len() < MAX_CACHED_PAIRS {
                    self.distances.insert(key, Some(0));
                }
                return Some(0);
            }
        }

        self.calculations += 1;
        let distance =
            self.scratch
                .distance_with_limit(first.as_bytes(), second.as_bytes(), max_distance);
        if let (Some(first_id), Some(second_id)) = (first_id, second_id) {
            if self.distances.len() < MAX_CACHED_PAIRS {
                self.distances
                    .insert((first_id, second_id, max_distance), distance);
            }
        }
        distance
    }
}

fn encode_qgrams(sequence: &[u8]) -> Vec<u32> {
    let mut qgrams: Vec<u32> = sequence
        .windows(QGRAM_LENGTH)
        .map(|gram| {
            gram.iter()
                .fold(0u32, |encoded, &byte| (encoded << 8) | u32::from(byte))
        })
        .collect();
    qgrams.sort_unstable();
    qgrams
}

fn qgram_excess(query: &[u32], target: &[u32]) -> usize {
    let (mut query_index, mut target_index, mut excess) = (0, 0, 0);
    while query_index < query.len() {
        let query_gram = query[query_index];
        let query_end = query[query_index..]
            .iter()
            .position(|&gram| gram != query_gram)
            .map_or(query.len(), |offset| query_index + offset);
        while target_index < target.len() && target[target_index] < query_gram {
            target_index += 1;
        }
        let target_count = if target_index < target.len() && target[target_index] == query_gram {
            target[target_index..]
                .iter()
                .take_while(|&&gram| gram == query_gram)
                .count()
        } else {
            0
        };
        excess += (query_end - query_index).saturating_sub(target_count);
        query_index = query_end;
    }
    excess
}

/// Cluster nodes by their exact HW sequence distance without building a distance matrix.
pub fn cluster_nodes_by_sequence(
    graph: &Graph<GeneNode, GeneEdge>,
    nodes: &[NodeIndex],
    threshold: f64,
    dna: bool,
    cache: &mut PairwiseDistanceCache,
) -> Vec<Vec<NodeIndex>> {
    let sequences: Vec<&str> = nodes
        .iter()
        .map(|&node_idx| {
            let node = &graph[node_idx];
            let sequence = if dna {
                node.dna.get(node.max_len_id)
            } else {
                None
            }
            .or_else(|| node.protein.get(node.max_len_id));
            sequence.map(String::as_str).unwrap_or("")
        })
        .collect();

    let sequence_ids: Vec<Option<u32>> = sequences
        .iter()
        .map(|sequence| cache.sequence_id(sequence))
        .collect();

    single_linkage_clustering_by(nodes, |i, j| {
        let seq1 = sequences[i];
        let seq2 = sequences[j];
        if seq1.is_empty() || seq2.is_empty() {
            return 1.0 < threshold;
        }

        let max_length = seq1.len().max(seq2.len());
        let max_dist = (max_length as f64 * (1.0 - threshold)).ceil() as usize;
        let minimum_length_distance = seq1.len().saturating_sub(seq2.len());
        let distance = if minimum_length_distance > max_dist {
            1.0
        } else {
            let Some(edit_distance) =
                cache.distance(sequence_ids[i], sequence_ids[j], seq1, seq2, max_dist)
            else {
                return 1.0 < threshold;
            };
            edit_distance as f64 / max_length as f64
        };
        distance < threshold
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::{align_sequences, calculate_identity};
    use crate::single_linkage::single_linkage_clustering;

    #[test]
    fn test_pwdist_empty() {
        let graph: Graph<GeneNode, GeneEdge> = Graph::default();
        let nodes = vec![];
        let result = cluster_nodes_by_sequence(
            &graph,
            &nodes,
            0.95,
            false,
            &mut PairwiseDistanceCache::default(),
        );
        assert!(result.is_empty());
    }

    #[test]
    fn test_pwdist_reuses_cached_pair_distances() {
        let mut graph: Graph<GeneNode, GeneEdge> = Graph::default();
        let first = graph.add_node(GeneNode {
            protein: vec!["ACGT".to_string()],
            ..Default::default()
        });
        let second = graph.add_node(GeneNode {
            protein: vec!["ACGA".to_string()],
            ..Default::default()
        });
        let nodes = [first, second];
        let mut cache = PairwiseDistanceCache::default();

        let first_result = cluster_nodes_by_sequence(&graph, &nodes, 0.7, false, &mut cache);
        assert_eq!(cache.calculation_count(), 1);
        assert_eq!(cache.hit_count(), 0);

        let second_result = cluster_nodes_by_sequence(&graph, &nodes, 0.7, false, &mut cache);
        assert_eq!(cache.calculation_count(), 1);
        assert_eq!(cache.hit_count(), 1);
        assert_eq!(first_result, second_result);
    }

    #[test]
    fn test_lazy_clustering_matches_materialized_reference() {
        let sequences = ["ACGT", "ACGA", "TTTT", "ACG", "ACGT"];
        let mut graph: Graph<GeneNode, GeneEdge> = Graph::default();
        let nodes: Vec<_> = sequences
            .iter()
            .map(|sequence| {
                graph.add_node(GeneNode {
                    protein: vec![sequence.to_string()],
                    ..Default::default()
                })
            })
            .collect();

        for threshold in [0.2, 0.4, 0.7, 0.95] {
            let mut reference = vec![vec![0.0; nodes.len()]; nodes.len()];
            for i in 0..nodes.len() {
                for j in (i + 1)..nodes.len() {
                    let first = sequences[i].as_bytes();
                    let second = sequences[j].as_bytes();
                    let max_distance =
                        (first.len().max(second.len()) as f64 * (1.0 - threshold)).ceil() as usize;
                    let alignment = align_sequences(first, second, max_distance);
                    let distance = if alignment.edit_distance == usize::MAX {
                        1.0
                    } else {
                        let identity = calculate_identity(
                            alignment.edit_distance,
                            first.len().max(second.len()),
                        );
                        1.0 - identity
                    };
                    reference[i][j] = distance;
                    reference[j][i] = distance;
                }
            }

            let mut cache = PairwiseDistanceCache::default();
            let actual = cluster_nodes_by_sequence(&graph, &nodes, threshold, false, &mut cache);
            let expected = single_linkage_clustering(&reference, &nodes, threshold);
            assert_eq!(canonicalize(actual), canonicalize(expected));
        }
    }

    #[test]
    fn test_connected_pairs_are_not_recomputed() {
        let mut graph: Graph<GeneNode, GeneEdge> = Graph::default();
        let nodes: Vec<_> = (0..6)
            .map(|index| {
                let mut sequence = vec![b'A'; 20];
                if index > 0 {
                    sequence[index] = b'C';
                }
                graph.add_node(GeneNode {
                    protein: vec![String::from_utf8(sequence).unwrap()],
                    ..Default::default()
                })
            })
            .collect();
        let mut cache = PairwiseDistanceCache::default();

        let clusters = cluster_nodes_by_sequence(&graph, &nodes, 0.7, false, &mut cache);

        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0], nodes);
        assert_eq!(cache.calculation_count(), nodes.len() - 1);
    }

    #[test]
    fn test_qgram_filter_never_rejects_an_in_range_hw_match() {
        let mut sequences = Vec::new();
        for length in 1..=6 {
            for value in 0..(1usize << length) {
                let sequence: Vec<u8> = (0..length)
                    .map(|bit| if value & (1 << bit) == 0 { b'A' } else { b'C' })
                    .collect();
                sequences.push(sequence);
            }
        }

        for query in &sequences {
            for target in &sequences {
                let mut scratch = crate::alignment::HwEditScratch::default();
                let distance = scratch.distance(query, target);
                for limit in 0..=3 {
                    if distance <= limit {
                        assert!(
                            qgram_excess(&encode_qgrams(query), &encode_qgrams(target))
                                <= QGRAM_LENGTH * limit,
                            "false rejection for {:?} vs {:?} at {}",
                            query,
                            target,
                            limit
                        );
                    }
                }
            }
        }
    }

    fn canonicalize(mut clusters: Vec<Vec<NodeIndex>>) -> Vec<Vec<NodeIndex>> {
        for cluster in &mut clusters {
            cluster.sort();
        }
        clusters.sort_by_key(|cluster| cluster.first().copied());
        clusters
    }
}
