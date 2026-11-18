//! Global candidate preclustering and sparse centroid similarity.

use crate::cdhit::{parse_cdhit_clusters, run_cluster, CdhitConfig};
use crate::network::{GeneEdge, GeneNode, Graph};
use anyhow::{bail, Context, Result};
use bio::pattern_matching::myers::{long::Myers, MyersBuilder};
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;
use tempfile::TempDir;

pub struct FamilyIndex {
    pub centroid_to_index: HashMap<String, usize>,
    pub close: Vec<HashSet<usize>>,
}

impl FamilyIndex {
    pub fn build(
        graph: &Graph<GeneNode, GeneEdge>,
        outdir: &str,
        threshold: f64,
        length_ratio: f64,
        threads: usize,
        dna: bool,
    ) -> Result<Self> {
        if threads == 0 || !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
            bail!("Family clustering requires positive threads and an identity between 0 and 1");
        }
        let started = Instant::now();
        let mut names = Vec::new();
        let mut sequences = Vec::new();
        let mut centroid_to_index = HashMap::new();
        for node in graph.node_weights() {
            let source = if dna { &node.dna } else { &node.protein };
            if source.len() < node.centroid.len() {
                bail!(
                    "Missing {} sequences for graph centroids",
                    if dna { "DNA" } else { "protein" }
                );
            }
            for (centroid, sequence) in node.centroid.iter().zip(source) {
                if centroid.is_empty() || sequence.is_empty() {
                    bail!("Empty centroid or sequence encountered in the pangenome graph");
                }
                if let Some(&index) = centroid_to_index.get(centroid) {
                    sequences[index] = sequence.clone();
                } else {
                    centroid_to_index.insert(centroid.clone(), names.len());
                    names.push(centroid.clone());
                    sequences.push(sequence.clone());
                }
            }
        }
        if names.is_empty() {
            return Ok(Self {
                centroid_to_index,
                close: Vec::new(),
            });
        }

        let temp = TempDir::new_in(outdir)?;
        let input = temp.path().join("centroids.fasta");
        let mut writer = BufWriter::new(File::create(&input)?);
        for (name, sequence) in names.iter().zip(&sequences) {
            writeln!(writer, ">{}\n{}", name, sequence)?;
        }
        writer.flush()?;
        drop(writer);
        let thresholds: &[f64] = if dna {
            &[0.99, 0.98, 0.95, 0.9]
        } else {
            &[0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5]
        };
        let mut clusters: Vec<Vec<usize>> = (0..names.len()).map(|index| vec![index]).collect();
        let mut current_input = input;
        for (iteration, &identity) in thresholds.iter().enumerate() {
            let output = temp.path().join(format!("pass_{}.fasta", iteration));
            run_cluster(
                path_string(&current_input)?,
                path_string(&output)?,
                &CdhitConfig {
                    identity,
                    length_ratio,
                    threads,
                    dna,
                    accurate: !dna,
                    quiet: true,
                },
            )?;
            let parsed = parse_cdhit_clusters(&format!("{}.clstr", path_string(&output)?))?;
            let mut representative_to_cluster = HashMap::new();
            for (cluster, members) in parsed {
                for member in members {
                    representative_to_cluster.insert(member, cluster.clone());
                }
            }
            let mut combined: BTreeMap<String, Vec<usize>> = BTreeMap::new();
            for cluster in clusters {
                let representative_cluster = cluster
                    .iter()
                    .filter_map(|&index| representative_to_cluster.get(&names[index]))
                    .last()
                    .context("CD-HIT lost a cluster representative during global preclustering")?;
                combined
                    .entry(representative_cluster.clone())
                    .or_default()
                    .extend(cluster);
            }
            clusters = combined.into_values().collect();
            current_input = output;
        }

        for cluster in &mut clusters {
            cluster.sort_by_key(|&index| sequences[index].len());
        }
        let pair_count: usize = clusters
            .iter()
            .map(|cluster| cluster.len() * cluster.len().saturating_sub(1) / 2)
            .sum();
        let builder = equality_builder(dna);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()?;
        let accepted: Vec<Vec<(usize, usize)>> = pool.install(|| {
            clusters
                .par_iter()
                .map(|cluster| {
                    cluster
                        .par_iter()
                        .enumerate()
                        .map(|(offset, &first)| {
                            let mut edges = Vec::new();
                            if offset + 1 == cluster.len() {
                                return edges;
                            }
                            let query = PreparedQuery::new(&sequences[first], dna, &builder);
                            let limit = query.threshold_distance(threshold);
                            for &second in &cluster[offset + 1..] {
                                if threshold == 0.0
                                    || query.matches_with_limit(&sequences[second], limit)
                                {
                                    edges.push((first, second));
                                }
                            }
                            edges
                        })
                        .collect::<Vec<_>>()
                        .into_iter()
                        .flatten()
                        .collect()
                })
                .collect()
        });
        let mut close = vec![HashSet::new(); names.len()];
        let mut accepted_count = 0;
        for (first, second) in accepted.into_iter().flatten() {
            close[first].insert(second);
            close[second].insert(first);
            accepted_count += 1;
        }
        log::info!(
            "Global {} similarity: {} centroids, {} candidate pairs, {} accepted pairs, {} threads, {:.3}s",
            if dna { "DNA" } else { "protein" }, names.len(), pair_count, accepted_count, threads,
            started.elapsed().as_secs_f64()
        );
        Ok(Self {
            centroid_to_index,
            close,
        })
    }

    pub fn clusters(
        &self,
        graph: &Graph<GeneNode, GeneEdge>,
        nodes: &[NodeIndex],
    ) -> Result<Vec<Vec<NodeIndex>>> {
        let mut centroid_nodes: HashMap<usize, Vec<usize>> = HashMap::new();
        for (position, &node) in nodes.iter().enumerate() {
            for centroid in &graph[node].centroid {
                let index = *self.centroid_to_index.get(centroid).with_context(|| {
                    format!("Missing centroid in family similarity index: {}", centroid)
                })?;
                centroid_nodes.entry(index).or_default().push(position);
            }
        }
        let mut parent: Vec<_> = (0..nodes.len()).collect();
        for (&centroid, positions) in &centroid_nodes {
            for &position in &positions[1..] {
                join(&mut parent, positions[0], position);
            }
            if self.close[centroid].len() <= centroid_nodes.len() {
                for neighbor in &self.close[centroid] {
                    if let Some(other_positions) = centroid_nodes.get(neighbor) {
                        join(&mut parent, positions[0], other_positions[0]);
                    }
                }
            } else {
                for (neighbor, other_positions) in &centroid_nodes {
                    if self.close[centroid].contains(neighbor) {
                        join(&mut parent, positions[0], other_positions[0]);
                    }
                }
            }
        }
        let mut result = Vec::<Vec<NodeIndex>>::new();
        let mut root_cluster = HashMap::new();
        for (position, &node) in nodes.iter().enumerate() {
            let root = root(&mut parent, position);
            let next = result.len();
            let cluster = *root_cluster.entry(root).or_insert_with(|| {
                result.push(Vec::new());
                next
            });
            result[cluster].push(node);
        }
        Ok(result)
    }
}

fn root(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}

fn join(parent: &mut [usize], first: usize, second: usize) {
    let first = root(parent, first);
    let second = root(parent, second);
    parent[second] = first;
}

fn path_string(path: &Path) -> Result<&str> {
    path.to_str().context("CD-HIT requires a UTF-8 path")
}

pub fn reference_identity(first: &str, second: &str, dna: bool) -> f64 {
    let (query, target) = if first.len() > second.len() {
        (second, first)
    } else {
        (first, second)
    };
    if query.is_empty() || target.is_empty() {
        return 0.0;
    }
    let builder = equality_builder(dna);
    PreparedQuery::new(query, dna, &builder).identity(target)
}

struct PreparedQuery {
    forward: Myers<u64>,
    reverse: Option<Myers<u64>>,
    length: usize,
}

impl PreparedQuery {
    fn new(query: &str, dna: bool, builder: &MyersBuilder) -> Self {
        Self {
            forward: builder.build_long_64(query.as_bytes()),
            reverse: dna.then(|| {
                builder.build_long_64(crate::prokka::reverse_complement(query).as_bytes())
            }),
            length: query.len(),
        }
    }

    fn identity(&self, target: &str) -> f64 {
        let score = |matcher: &Myers<u64>| {
            matcher
                .find_all_end(target.as_bytes(), self.length / 2)
                .map(|(_, distance)| 1.0 - distance as f64 / self.length as f64)
                .fold(0.0, f64::max)
        };
        let forward = score(&self.forward);
        if forward == 1.0 {
            return forward;
        }
        self.reverse
            .as_ref()
            .map_or(forward, |reverse| forward.max(score(reverse)))
    }

    fn threshold_distance(&self, threshold: f64) -> usize {
        let passes = |distance: usize| 1.0 - distance as f64 / self.length as f64 >= threshold;
        let (mut low, mut high) = (0, self.length / 2 + 1);
        while low < high {
            let mid = low + (high - low) / 2;
            if passes(mid) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        low - 1
    }

    fn matches_with_limit(&self, target: &str, limit: usize) -> bool {
        self.forward
            .find_all_end(target.as_bytes(), limit)
            .next()
            .is_some()
            || self.reverse.as_ref().is_some_and(|reverse| {
                reverse
                    .find_all_end(target.as_bytes(), limit)
                    .next()
                    .is_some()
            })
    }
}

fn equality_builder(dna: bool) -> MyersBuilder {
    let pairs: Vec<(u8, u8)> = if dna {
        b"ACGT".iter().map(|&base| (base, b'N')).collect()
    } else {
        b"*ACBEDGFIHKMLNQPSRTWVYXZ"
            .iter()
            .map(|&base| (base, b'X'))
            .chain([(b'D', b'B'), (b'N', b'B'), (b'E', b'Z'), (b'Q', b'Z')])
            .collect()
    };
    let mut equalities: HashMap<u8, Vec<u8>> = HashMap::new();
    for (first, second) in pairs {
        equalities.entry(first).or_default().push(second);
        equalities.entry(second).or_default().push(first);
    }
    let mut builder = MyersBuilder::new();
    for (base, equivalents) in equalities {
        builder.ambig(base, equivalents);
    }
    builder
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_uses_shorter_sequence_and_ambiguities() {
        assert_eq!(reference_identity("ABC", "ZZABCZZ", false), 1.0);
        assert_eq!(reference_identity("DNEQ", "BBZZ", false), 1.0);
        assert_eq!(reference_identity("AAAA", "XXXX", false), 1.0);
        assert_eq!(reference_identity("ACGT", "NNNN", true), 1.0);
        assert_eq!(reference_identity("AAAA", "TTTT", true), 1.0);
        assert_eq!(reference_identity("AAAA", "TTTT", false), 0.0);
    }

    #[test]
    fn prepared_query_reuse_preserves_independent_target_scores() {
        for dna in [false, true] {
            let builder = equality_builder(dna);
            let query = PreparedQuery::new("ACGT", dna, &builder);
            for target in ["ACGT", "TTTTTT", "ACGTTT", "XXXX", "NNNN", "ACGT"] {
                assert_eq!(
                    query.identity(target),
                    reference_identity("ACGT", target, dna)
                );
            }
        }
    }

    #[test]
    fn bounded_threshold_matches_full_identity_at_boundaries() {
        for dna in [false, true] {
            let builder = equality_builder(dna);
            for query in ["A", "ACGT", "DNEQ", "ACGTACGT", "NNNN"] {
                let prepared = PreparedQuery::new(query, dna, &builder);
                for target in ["A", "ACGT", "TTTT", "BBZZ", "ACGTACGA", "NNNN", "XXXX", ""] {
                    let identity = prepared.identity(target);
                    for threshold in [
                        0.0,
                        0.25,
                        0.5,
                        0.75,
                        0.9,
                        1.0,
                        identity,
                        (identity + f64::EPSILON).min(1.0),
                    ] {
                        assert_eq!(
                            threshold == 0.0
                                || prepared.matches_with_limit(
                                    target,
                                    prepared.threshold_distance(threshold)
                                ),
                            identity >= threshold,
                            "dna={dna} query={query} target={target} threshold={threshold}"
                        );
                    }
                    for dna in [false, true] {
                        let builder = equality_builder(dna);
                        for length in [63, 64, 65, 127, 128, 129] {
                            let query: String = (0..length)
                                .map(|i| "ACGT".as_bytes()[(i * 3 + i / 7) % 4] as char)
                                .collect();
                            let prepared = PreparedQuery::new(&query, dna, &builder);
                            let mut target = query.clone().into_bytes();
                            for position in (3..length).step_by(9) {
                                target[position] = b'N';
                            }
                            let target = String::from_utf8(target).unwrap();
                            let identity = prepared.identity(&target);
                            for threshold in [
                                0.0,
                                0.5,
                                0.9,
                                1.0,
                                identity,
                                (identity + f64::EPSILON).min(1.0),
                            ] {
                                assert_eq!(
                                    threshold == 0.0
                                        || prepared.matches_with_limit(
                                            &target,
                                            prepared.threshold_distance(threshold)
                                        ),
                                    identity >= threshold,
                                    "dna={dna} length={length} threshold={threshold}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn sparse_clustering_connects_all_centroids_of_merged_nodes() {
        let mut graph = Graph::default();
        let first = graph.add_node(GeneNode {
            centroid: vec!["a".into(), "b".into()],
            ..Default::default()
        });
        let second = graph.add_node(GeneNode {
            centroid: vec!["c".into()],
            ..Default::default()
        });
        let third = graph.add_node(GeneNode {
            centroid: vec!["d".into()],
            ..Default::default()
        });
        let index = FamilyIndex {
            centroid_to_index: ["a", "b", "c", "d"]
                .iter()
                .enumerate()
                .map(|(index, name)| (name.to_string(), index))
                .collect(),
            close: vec![
                HashSet::new(),
                HashSet::from([2]),
                HashSet::from([1]),
                HashSet::new(),
            ],
        };
        assert_eq!(
            index.clusters(&graph, &[first, second, third]).unwrap(),
            vec![vec![first, second], vec![third]]
        );
    }

    #[test]
    fn clustering_ignores_nonlocal_neighbours_and_keeps_node_order() {
        let mut graph = Graph::default();
        let a = graph.add_node(GeneNode {
            centroid: vec!["a".into()],
            ..Default::default()
        });
        let b = graph.add_node(GeneNode {
            centroid: vec!["b".into()],
            ..Default::default()
        });
        let c = graph.add_node(GeneNode {
            centroid: vec!["c".into()],
            ..Default::default()
        });
        let mut close = vec![HashSet::new(); 20];
        close[0].extend(1..20);
        close[1].insert(0);
        let mut index = FamilyIndex {
            centroid_to_index: HashMap::from([("a".into(), 0), ("b".into(), 1), ("c".into(), 2)]),
            close,
        };
        assert_eq!(
            index.clusters(&graph, &[c, a, b]).unwrap(),
            vec![vec![c, a, b]]
        );
        index.close[0].remove(&2);
        assert_eq!(
            index.clusters(&graph, &[c, a, b]).unwrap(),
            vec![vec![c], vec![a, b]]
        );
    }
}
