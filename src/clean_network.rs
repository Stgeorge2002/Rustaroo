//! Network cleaning algorithms - family collapse and graph cleanup
//!
//! Port of core functions from rustaroo/clean_network.py

use crate::family_index::FamilyIndex;
use crate::network::Graph;
use crate::network::{GeneEdge, GeneNode};
use anyhow::Result;
use petgraph::graph::NodeIndex;
use petgraph::visit::{EdgeRef, NodeIndexable};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::time::Instant;

pub fn clean_misassembly_edges(graph: &mut Graph<GeneNode, GeneEdge>, threshold: usize) {
    let bad: Vec<_> = graph
        .edge_indices()
        .filter(|&edge| {
            let (first, second) = graph
                .edge_endpoints(edge)
                .expect("Existing edge has endpoints");
            let support = graph[edge].size;
            support < threshold
                && (graph[first].has_end
                    || graph[second].has_end
                    || (support as f64) < 0.05 * graph[first].size.min(graph[second].size) as f64)
        })
        .collect();
    for edge in bad {
        graph.remove_edge(edge);
    }
}

/// Trim low support trailing ends (genes at contig ends with low support)
pub fn trim_low_support_trailing_ends(
    graph: &mut Graph<GeneNode, GeneEdge>,
    min_support: usize,
    max_recursive: usize,
) -> Result<()> {
    let mut candidates: Vec<_> = graph.node_indices().collect();
    for _iteration in 0..max_recursive {
        let bad_nodes: Vec<_> = candidates
            .into_iter()
            .filter(|&node| {
                graph.node_weight(node).is_some()
                    && graph[node].size < min_support
                    && graph.edges(node).count() <= 1
            })
            .collect();
        if bad_nodes.is_empty() {
            break;
        }
        let mut affected = Vec::new();
        for node in bad_nodes {
            affected.extend(graph.neighbors(node));
            graph.remove_node(node);
        }
        affected.sort_unstable();
        affected.dedup();
        candidates = affected;
    }

    Ok(())
}

/// Modified BFS to explore neighborhood up to depth limit
pub fn mod_bfs_edges(
    graph: &Graph<GeneNode, GeneEdge>,
    source: NodeIndex,
    depth_limit: Option<usize>,
) -> Vec<(NodeIndex, NodeIndex, usize)> {
    let mut scratch = BfsScratch::default();
    mod_bfs_edges_with_scratch(graph, source, depth_limit, true, &mut scratch)
}

#[derive(Default)]
struct BfsScratch {
    visited_generation: Vec<u32>,
    generation: u32,
    queue: VecDeque<(NodeIndex, usize)>,
    edges: Vec<(NodeIndex, NodeIndex, usize)>,
    reached: Vec<NodeIndex>,
}

fn mod_bfs_edges_with_scratch(
    graph: &Graph<GeneNode, GeneEdge>,
    source: NodeIndex,
    depth_limit: Option<usize>,
    record_edges: bool,
    scratch: &mut BfsScratch,
) -> Vec<(NodeIndex, NodeIndex, usize)> {
    scratch.queue.clear();
    scratch.edges.clear();
    scratch.reached.clear();
    if scratch.visited_generation.len() < graph.node_bound() {
        scratch.visited_generation.resize(graph.node_bound(), 0);
    }
    scratch.generation = scratch.generation.wrapping_add(1);
    if scratch.generation == 0 {
        scratch.visited_generation.fill(0);
        scratch.generation = 1;
    }
    scratch.visited_generation[source.index()] = scratch.generation;
    let depth_limit = depth_limit.unwrap_or(graph.node_count());
    scratch.queue.push_back((source, depth_limit));

    while let Some((parent, depth_now)) = scratch.queue.pop_front() {
        for child in graph.neighbors(parent) {
            if scratch.visited_generation[child.index()] != scratch.generation {
                scratch.visited_generation[child.index()] = scratch.generation;
                if record_edges {
                    scratch.edges.push((parent, child, depth_now));
                }
                scratch.reached.push(child);
                if depth_now > 1 {
                    scratch.queue.push_back((child, depth_now - 1));
                }
            }
        }
    }

    scratch.edges.clone()
}

#[cfg(test)]
mod bfs_tests {
    use super::*;

    #[test]
    fn misassembly_cleaning_respects_support_and_fraction_boundaries() {
        let mut graph = Graph::<GeneNode, GeneEdge>::default();
        let large = graph.add_node(GeneNode {
            size: 100,
            ..Default::default()
        });
        let end = graph.add_node(GeneNode {
            size: 100,
            has_end: true,
            ..Default::default()
        });
        let equal_fraction = graph.add_node(GeneNode {
            size: 100,
            ..Default::default()
        });
        let low_fraction = graph.add_node(GeneNode {
            size: 100,
            ..Default::default()
        });
        let at_threshold = graph.add_node(GeneNode {
            size: 100,
            has_end: true,
            ..Default::default()
        });
        for (node, size) in [
            (end, 5),
            (equal_fraction, 5),
            (low_fraction, 4),
            (at_threshold, 6),
        ] {
            graph.add_edge(
                large,
                node,
                GeneEdge {
                    size,
                    members: BTreeSet::new(),
                },
            );
        }
        clean_misassembly_edges(&mut graph, 6);
        assert!(graph.find_edge(large, end).is_none());
        assert!(graph.find_edge(large, equal_fraction).is_some());
        assert!(graph.find_edge(large, low_fraction).is_none());
        assert!(graph.find_edge(large, at_threshold).is_some());
    }

    #[test]
    fn reusable_bfs_preserves_depth_limited_reachability() {
        let mut graph = Graph::<GeneNode, GeneEdge>::default();
        let nodes: Vec<_> = (0..5)
            .map(|_| graph.add_node(GeneNode::default()))
            .collect();
        for (first, second) in [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4)] {
            graph.add_edge(
                nodes[first],
                nodes[second],
                GeneEdge {
                    size: 1,
                    members: BTreeSet::new(),
                },
            );
        }

        let mut scratch = BfsScratch::default();
        let first = mod_bfs_edges_with_scratch(&graph, nodes[0], Some(2), true, &mut scratch);
        assert_eq!(first.len(), 3);
        assert!(first.contains(&(nodes[0], nodes[1], 2)));
        assert!(first.contains(&(nodes[0], nodes[2], 2)));
        assert!(
            first.contains(&(nodes[1], nodes[3], 1)) || first.contains(&(nodes[2], nodes[3], 1))
        );

        let second = mod_bfs_edges_with_scratch(&graph, nodes[3], Some(1), true, &mut scratch);
        assert_eq!(second.len(), 3);
    }
}

/// Cluster neighbors by sequence similarity (simplified single linkage)
pub fn cluster_neighbors(
    graph: &Graph<GeneNode, GeneEdge>,
    neighbours: &[NodeIndex],
    similarity_threshold: f64,
) -> Vec<Vec<NodeIndex>> {
    let mut clusters: Vec<Vec<NodeIndex>> = Vec::new();

    for &node in neighbours {
        let mut added = false;

        // Try to add to existing cluster
        for cluster in &mut clusters {
            let first_node = cluster[0];

            // Check if similar enough (simplified - in reality uses alignment)
            if are_nodes_similar(graph, node, first_node, similarity_threshold) {
                cluster.push(node);
                added = true;
                break;
            }
        }

        if !added {
            clusters.push(vec![node]);
        }
    }

    clusters
}

/// Check if two nodes are similar (simplified version)
fn are_nodes_similar(
    graph: &Graph<GeneNode, GeneEdge>,
    node1: NodeIndex,
    node2: NodeIndex,
    threshold: f64,
) -> bool {
    let seq1 = &graph[node1].protein[0];
    let seq2 = &graph[node2].protein[0];

    // Simple length-based similarity check (real version uses alignment)
    let len_diff = (seq1.len() as f64 - seq2.len() as f64).abs();
    let max_len = seq1.len().max(seq2.len()) as f64;

    if max_len == 0.0 {
        return false;
    }

    let similarity = 1.0 - (len_diff / max_len);
    similarity >= threshold
}

/// Collapse gene families - merge similar neighboring gene clusters
pub fn collapse_families(
    graph: &mut Graph<GeneNode, GeneEdge>,
    outdir: &str,
    family_threshold: f64,
    len_dif_percent: f64,
    depths: &[usize],
    n_cpu: usize,
    quiet: bool,
) -> Result<()> {
    let seqid_to_centroid = graph
        .node_weights()
        .flat_map(|node| {
            node.seq_ids
                .iter()
                .map(|seq_id| (seq_id.clone(), node.centroid[0].clone()))
        })
        .collect();
    let index = FamilyIndex::build(
        graph,
        outdir,
        family_threshold,
        len_dif_percent,
        n_cpu,
        false,
    )?;
    collapse_families_indexed(
        graph,
        &seqid_to_centroid,
        &index,
        depths,
        quiet,
        false,
        0.01,
    )
}

pub fn collapse_families_indexed(
    graph: &mut Graph<GeneNode, GeneEdge>,
    seqid_to_centroid: &HashMap<String, String>,
    index: &FamilyIndex,
    depths: &[usize],
    quiet: bool,
    dna: bool,
    outlier_support: f64,
) -> Result<()> {
    let mut node_members: HashMap<NodeIndex, HashMap<usize, HashSet<usize>>> = HashMap::new();
    for node in graph.node_indices() {
        let mut members = HashMap::<usize, HashSet<usize>>::new();
        for seq_id in &graph[node].seq_ids {
            let member = seq_id.split('_').next().unwrap_or("").parse::<usize>()?;
            let centroid = if seq_id.contains("refound") {
                &graph[node].long_centroid_id.1
            } else {
                seqid_to_centroid
                    .get(seq_id)
                    .ok_or_else(|| anyhow::anyhow!("Missing centroid for sequence {}", seq_id))?
            };
            let centroid_index = index
                .centroid_to_index
                .get(centroid)
                .ok_or_else(|| anyhow::anyhow!("Missing indexed centroid {}", centroid))?;
            members.entry(member).or_default().insert(*centroid_index);
        }
        node_members.insert(node, members);
    }
    let mut bfs_scratch = BfsScratch::default();

    for &depth in depths {
        let depth_started = Instant::now();
        if !quiet {
            println!("Processing depth: {}", depth);
        }

        let mut search_space: BTreeSet<_> = graph.node_indices().collect();
        let mut iteration_num = 1;

        while !search_space.is_empty() {
            if !quiet {
                println!("Iteration: {}", iteration_num);
            }
            iteration_num += 1;
            let iteration_started = Instant::now();
            let mut centers_examined = 0usize;
            let mut neighborhood_nodes = 0usize;
            let mut merges = 0usize;

            let temp_node_list: Vec<_> = search_space.iter().cloned().collect();
            let mut removed_nodes = HashSet::new();

            for node in temp_node_list {
                if removed_nodes.contains(&node) || graph.node_weight(node).is_none() {
                    continue;
                }

                let degree: usize = graph
                    .edges(node)
                    .take(3)
                    .map(|edge| if edge.source() == edge.target() { 2 } else { 1 })
                    .sum();
                if degree <= 2 {
                    if !removed_nodes.contains(&node) {
                        search_space.remove(&node);
                    }
                    removed_nodes.insert(node);
                    continue;
                }
                centers_examined += 1;

                // Find neighboring nodes
                mod_bfs_edges_with_scratch(graph, node, Some(depth), false, &mut bfs_scratch);
                let mut neighbours = std::mem::take(&mut bfs_scratch.reached);
                neighbours.push(node);
                neighborhood_nodes += neighbours.len();

                let clusters = index.clusters(graph, &neighbours)?;
                bfs_scratch.reached = neighbours;

                for cluster in clusters {
                    if cluster.len() <= 1 {
                        continue;
                    }

                    // Check for conflicts (same genome in multiple nodes)
                    let mut seen = HashSet::new();
                    let mut no_conflict = true;

                    for &n in &cluster {
                        if graph.node_weight(n).is_none() {
                            no_conflict = false;
                            break;
                        }
                        for &member in &graph[n].members {
                            if seen.contains(&member) {
                                no_conflict = false;
                                break;
                            }
                            seen.insert(member);
                        }
                        if !no_conflict {
                            break;
                        }
                    }

                    if no_conflict {
                        merges += 1;
                        merge_family_cluster(
                            graph,
                            &cluster,
                            &mut node_members,
                            &mut search_space,
                            &mut removed_nodes,
                            !dna,
                        )?;
                    } else {
                        let mut cluster = cluster;
                        cluster.sort_by_key(|&node| std::cmp::Reverse(graph[node].size));
                        let mut member_counts = HashMap::<usize, usize>::new();
                        for &node in &cluster {
                            for &member in &graph[node].members {
                                *member_counts.entry(member).or_default() += 1;
                            }
                        }
                        let singleton_fraction =
                            member_counts.values().filter(|&&count| count == 1).count() as f64
                                / member_counts.len().max(1) as f64;
                        while !cluster.is_empty() {
                            let first = cluster[0];
                            let mut subset = vec![first];
                            for &second in &cluster[1..] {
                                if graph[first].members.is_disjoint(&graph[second].members)
                                    || (singleton_fraction >= outlier_support
                                        && split_gene_merge_allowed(
                                            graph,
                                            index,
                                            &node_members,
                                            first,
                                            second,
                                        ))
                                {
                                    subset.push(second);
                                }
                            }
                            if subset.len() > 1 {
                                for component in index.clusters(graph, &subset)? {
                                    if component.len() > 1 {
                                        merges += 1;
                                        merge_family_cluster(
                                            graph,
                                            &component,
                                            &mut node_members,
                                            &mut search_space,
                                            &mut removed_nodes,
                                            !dna,
                                        )?;
                                    }
                                }
                            }
                            let subset: HashSet<_> = subset.into_iter().collect();
                            cluster.retain(|node| !subset.contains(node));
                        }
                    }
                }
                if !removed_nodes.contains(&node) {
                    search_space.remove(&node);
                }
            }

            log::debug!(
                "Cleaning depth {} iteration {}: {} centers, {} neighborhood nodes, {} merges, {:.3}s",
                depth,
                iteration_num - 1,
                centers_examined,
                neighborhood_nodes,
                merges,
                iteration_started.elapsed().as_secs_f64()
            );

            if removed_nodes.is_empty() {
                break;
            }
        }

        log::debug!(
            "Cleaning depth {} completed in {:.3}s",
            depth,
            depth_started.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

fn merge_family_cluster(
    graph: &mut Graph<GeneNode, GeneEdge>,
    cluster: &[NodeIndex],
    node_members: &mut HashMap<NodeIndex, HashMap<usize, HashSet<usize>>>,
    search_space: &mut BTreeSet<NodeIndex>,
    removed_nodes: &mut HashSet<NodeIndex>,
    multi_centroid: bool,
) -> Result<()> {
    let mut members = HashMap::<usize, HashSet<usize>>::new();
    for node in cluster {
        if let Some(mapping) = node_members.remove(node) {
            for (member, centroids) in mapping {
                members.entry(member).or_default().extend(centroids);
            }
        }
        search_space.remove(node);
        removed_nodes.insert(*node);
    }
    let new_node =
        crate::merge_nodes::merge_node_cluster(graph, cluster, NodeIndex::new(0), multi_centroid)?;
    node_members.insert(new_node, members);
    search_space.insert(new_node);
    // Reused slots must not be processed as the old node in this iteration's snapshot.
    removed_nodes.insert(new_node);
    Ok(())
}

fn split_gene_merge_allowed(
    graph: &Graph<GeneNode, GeneEdge>,
    index: &FamilyIndex,
    node_members: &HashMap<NodeIndex, HashMap<usize, HashSet<usize>>>,
    first: NodeIndex,
    second: NodeIndex,
) -> bool {
    if graph[first]
        .centroid
        .iter()
        .any(|centroid| graph[second].centroid.contains(centroid))
    {
        return false;
    }
    let mut edges_seen = HashSet::new();
    let mut edge_member_count = HashMap::<usize, usize>::new();
    for node in [first, second] {
        for edge in graph.edges(node) {
            if !edges_seen.insert(edge.id()) {
                continue;
            }
            for &member in &edge.weight().members {
                let count = edge_member_count.entry(member).or_default();
                *count += 1;
                if *count > 3 {
                    return false;
                }
            }
        }
    }
    for member in graph[first].members.intersection(&graph[second].members) {
        for first_centroid in &node_members[&first][member] {
            for second_centroid in &node_members[&second][member] {
                if index.close[*first_centroid].contains(second_centroid) {
                    return false;
                }
            }
        }
    }
    true
}

/// Collapse paralogous gene copies
pub fn collapse_paralogs(
    graph: &mut Graph<GeneNode, GeneEdge>,
    centroid_context: &HashMap<String, Vec<(NodeIndex, usize)>>,
) -> Result<()> {
    let mut ordered_contexts: Vec<_> = centroid_context.iter().collect();
    ordered_contexts.sort_by_key(|(centroid, _)| *centroid);
    for (_, contexts) in ordered_contexts {
        let mut contexts = contexts.clone();
        contexts.sort();
        let mut member_paralogs = BTreeMap::<usize, Vec<NodeIndex>>::new();
        for &(node, member) in &contexts {
            member_paralogs.entry(member).or_default().push(node);
        }
        let Some((&reference_member, references)) =
            member_paralogs.iter().max_by(|first, second| {
                first
                    .1
                    .len()
                    .cmp(&second.1.len())
                    .then_with(|| second.0.cmp(first.0))
            })
        else {
            continue;
        };
        let mut clusters: Vec<Vec<NodeIndex>> = references.iter().map(|&node| vec![node]).collect();
        let mut members: Vec<_> = references
            .iter()
            .map(|_| BTreeSet::from([reference_member]))
            .collect();
        let targets: BTreeSet<_> = contexts
            .iter()
            .filter_map(|&(node, member)| (member != reference_member).then_some(node))
            .collect();
        let distances: Vec<_> = references
            .iter()
            .map(|&reference| graph_distances(graph, reference, &targets))
            .collect();
        let reference_contexts: Vec<_> = references
            .iter()
            .map(|&reference| weighted_context(graph, reference))
            .collect();
        for &(node, member) in &contexts {
            if member == reference_member {
                continue;
            }
            let mut best = None;
            let mut best_distance = usize::MAX;
            for (cluster, distance) in distances.iter().enumerate() {
                if !members[cluster].contains(&member) {
                    if let Some(&distance) = distance.get(&node) {
                        if distance < best_distance {
                            best_distance = distance;
                            best = Some(cluster);
                        }
                    }
                }
            }
            if best.is_none() {
                let context = weighted_context(graph, node);
                let mut best_score = f64::NEG_INFINITY;
                for (cluster, reference) in reference_contexts.iter().enumerate() {
                    if members[cluster].contains(&member) {
                        continue;
                    }
                    let score: f64 = context
                        .iter()
                        .filter_map(|(centroid, &depth)| {
                            reference.get(centroid).map(|&reference_depth| {
                                1.0 / (1.0 + depth.abs_diff(reference_depth) as f64)
                            })
                        })
                        .sum();
                    if score > best_score {
                        best_score = score;
                        best = Some(cluster);
                    }
                }
            }
            let best = best.ok_or_else(|| {
                anyhow::anyhow!("No non-conflicting paralog context for genome {}", member)
            })?;
            clusters[best].push(node);
            members[best].insert(member);
        }
        for cluster in clusters {
            if cluster.len() > 1 {
                crate::merge_nodes::merge_node_cluster(graph, &cluster, NodeIndex::new(0), true)?;
            }
        }
    }
    Ok(())
}

fn graph_distances(
    graph: &Graph<GeneNode, GeneEdge>,
    source: NodeIndex,
    targets: &BTreeSet<NodeIndex>,
) -> HashMap<NodeIndex, usize> {
    let mut distances = HashMap::from([(source, 0usize)]);
    let mut queue = VecDeque::from([source]);
    let mut remaining = targets.len() - usize::from(targets.contains(&source));
    while let Some(node) = queue.pop_front() {
        if remaining == 0 {
            break;
        }
        let distance = distances[&node] + 1;
        for neighbor in graph.neighbors(node) {
            if let std::collections::hash_map::Entry::Vacant(entry) = distances.entry(neighbor) {
                entry.insert(distance);
                queue.push_back(neighbor);
                if targets.contains(&neighbor) {
                    remaining -= 1;
                }
            }
        }
    }
    distances
}

fn weighted_context(graph: &Graph<GeneNode, GeneEdge>, node: NodeIndex) -> HashMap<String, usize> {
    let mut context = HashMap::new();
    for (_, neighbor, depth) in mod_bfs_edges(graph, node, Some(5)) {
        context.insert(graph[neighbor].centroid[0].clone(), depth);
    }
    context
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frontier_trim_matches_full_scans_for_rounds_and_parallel_edges() {
        let mut state = 17u64;
        for _ in 0..200 {
            let mut graph = Graph::default();
            let nodes: Vec<_> = (0..8)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    graph.add_node(GeneNode {
                        size: (state >> 32) as usize % 5,
                        ..Default::default()
                    })
                })
                .collect();
            for _ in 0..20 {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let a = nodes[(state >> 32) as usize % nodes.len()];
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let b = nodes[(state >> 32) as usize % nodes.len()];
                graph.add_edge(
                    a,
                    b,
                    GeneEdge {
                        size: 1,
                        members: BTreeSet::new(),
                    },
                );
            }
            for rounds in 0..=4 {
                let mut expected = graph.clone();
                for _ in 0..rounds {
                    let bad: Vec<_> = expected
                        .node_indices()
                        .filter(|&node| {
                            expected[node].size < 3 && expected.edges(node).count() <= 1
                        })
                        .collect();
                    if bad.is_empty() {
                        break;
                    }
                    for node in bad {
                        expected.remove_node(node);
                    }
                }
                let mut actual = graph.clone();
                trim_low_support_trailing_ends(&mut actual, 3, rounds).unwrap();
                assert_eq!(
                    actual.node_indices().collect::<Vec<_>>(),
                    expected.node_indices().collect::<Vec<_>>()
                );
                assert_eq!(
                    actual
                        .edge_indices()
                        .map(|edge| (edge, actual.edge_endpoints(edge)))
                        .collect::<Vec<_>>(),
                    expected
                        .edge_indices()
                        .map(|edge| (edge, expected.edge_endpoints(edge)))
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    #[test]
    fn limited_distances_match_full_bfs_for_reachable_and_unreachable_targets() {
        let mut graph = Graph::default();
        let nodes: Vec<_> = (0..8)
            .map(|_| graph.add_node(GeneNode::default()))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 3), (2, 4), (4, 5), (6, 7)] {
            graph.add_edge(
                nodes[a],
                nodes[b],
                GeneEdge {
                    size: 1,
                    members: BTreeSet::new(),
                },
            );
        }
        let all: BTreeSet<_> = nodes.iter().copied().collect();
        for source in &nodes {
            let full = graph_distances(&graph, *source, &all);
            for targets in [
                BTreeSet::new(),
                BTreeSet::from([nodes[2], nodes[3]]),
                BTreeSet::from([nodes[0], nodes[7]]),
            ] {
                let limited = graph_distances(&graph, *source, &targets);
                for target in targets {
                    assert_eq!(limited.get(&target), full.get(&target));
                }
            }
        }
    }

    #[test]
    fn test_trim_trailing_ends() {
        let mut graph = Graph::default();
        let node1 = graph.add_node(GeneNode {
            size: 2,
            ..Default::default()
        });
        let node2 = graph.add_node(GeneNode {
            size: 5,
            ..Default::default()
        });

        graph.add_edge(
            node1,
            node2,
            GeneEdge {
                size: 1,
                members: BTreeSet::new(),
            },
        );

        let initial_count = graph.node_count();
        trim_low_support_trailing_ends(&mut graph, 3, 2).unwrap();

        // Node with size 2 should be removed (below min_support of 3)
        assert!(graph.node_count() < initial_count);
    }
}
