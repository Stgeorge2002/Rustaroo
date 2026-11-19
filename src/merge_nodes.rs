//! Node merging utilities for graph manipulation
//!

use crate::network::{add_ordered_edge, ordered_edges, GeneEdge, GeneNode, Graph};
use anyhow::{Context, Result};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use std::collections::HashSet;

pub(crate) fn next_creation_order(graph: &Graph<GeneNode, GeneEdge>) -> usize {
    graph
        .node_weights()
        .map(|node| node.creation_order)
        .max()
        .map_or(0, |id| id + 1)
}

pub(crate) fn merge_node_in_creation_order(
    graph: &mut Graph<GeneNode, GeneEdge>,
    cluster: &[NodeIndex],
    next_order: &mut usize,
    multi_centroid: bool,
) -> Result<NodeIndex> {
    let node = merge_node_cluster(graph, cluster, NodeIndex::new(0), multi_centroid)?;
    graph[node].creation_order = *next_order;
    *next_order += 1;
    Ok(node)
}

/// Merge multiple nodes into a single new node
pub fn merge_node_cluster(
    graph: &mut Graph<GeneNode, GeneEdge>,
    cluster: &[NodeIndex],
    _new_node_id: NodeIndex,
    multi_centroid: bool,
) -> Result<NodeIndex> {
    if cluster.is_empty() {
        anyhow::bail!("Cannot merge empty cluster");
    }

    let mut cluster = cluster.to_vec();
    for &node in &cluster {
        if graph.node_weight(node).is_none() {
            anyhow::bail!("Cannot merge deleted graph node {}", node.index());
        }
    }
    cluster.sort_by_key(|&node| graph[node].size);
    let cluster_set: HashSet<_> = cluster.iter().copied().collect();
    let mut merged = GeneNode::default();
    let mut centroids = HashSet::new();
    let mut extra_dna = Vec::new();
    let mut extra_protein = Vec::new();

    for &node_idx in &cluster {
        let node = &graph[node_idx];
        merged.size += node.size;
        merged.members.extend(node.members.iter());
        merged.seq_ids.extend(node.seq_ids.iter().cloned());
        merged.lengths.extend(node.lengths.iter());
        for (index, centroid) in node.centroid.iter().enumerate() {
            let dna = node
                .dna
                .get(index)
                .with_context(|| format!("Missing DNA sequence for centroid {}", centroid))?;
            let protein = node
                .protein
                .get(index)
                .with_context(|| format!("Missing protein sequence for centroid {}", centroid))?;
            if centroids.insert(centroid.clone()) {
                merged.centroid.push(centroid.clone());
                merged.dna.push(dna.clone());
                merged.protein.push(protein.clone());
            }
        }
        extra_dna.extend(node.dna.iter().skip(node.centroid.len()).cloned());
        extra_protein.extend(node.protein.iter().skip(node.centroid.len()).cloned());

        if node.has_end {
            merged.has_end = true;
        }
        if node.paralog {
            merged.paralog = true;
        }

        merged.annotation.push_str(&format!(";{}", node.annotation));
        merged
            .description
            .push_str(&format!(";{}", node.description));
        merged.long_centroid_id = merged.long_centroid_id.max(node.long_centroid_id.clone());
        merged.merged_dna |= node.merged_dna || !multi_centroid;
    }

    append_unique(&mut merged.dna, extra_dna);
    append_unique(&mut merged.protein, extra_protein);
    merged.annotation = join_annotations(&merged.annotation);
    merged.description = join_annotations(&merged.description);
    merged.size = merged.members.len();
    if !merged.dna.is_empty() {
        merged.max_len_id = merged
            .dna
            .iter()
            .enumerate()
            .max_by_key(|(_, sequence)| sequence.len())
            .map(|(idx, _)| idx)
            .unwrap_or(0);
    }

    // Collect edges to neighbors (not within cluster)
    let mut neighbor_positions = std::collections::HashMap::<NodeIndex, usize>::new();
    let mut neighbor_edges = Vec::<(NodeIndex, GeneEdge)>::new();
    let mut adjacency = Vec::new();

    for &node_idx in &cluster {
        ordered_edges(graph, node_idx, &mut adjacency);
        for &(_, edge, neighbor) in &adjacency {
            if cluster_set.contains(&neighbor) {
                continue;
            }
            match neighbor_positions.entry(neighbor) {
                std::collections::hash_map::Entry::Vacant(slot) => {
                    slot.insert(neighbor_edges.len());
                    neighbor_edges.push((neighbor, graph[edge].clone()));
                }
                std::collections::hash_map::Entry::Occupied(slot) => {
                    let data = &mut neighbor_edges[*slot.get()].1;
                    data.members.extend(graph[edge].members.iter());
                    data.size = data.members.len();
                }
            }
        }
    }

    let actual_idx = graph.add_node(merged);

    // Add edges to neighbors
    for (neighbor, edge_data) in neighbor_edges {
        if graph.node_weight(neighbor).is_some() {
            add_ordered_edge(graph, actual_idx, neighbor, edge_data);
        }
    }
    for &node_idx in &cluster {
        graph.remove_node(node_idx);
    }

    fn append_unique(values: &mut Vec<String>, extra: Vec<String>) {
        let mut seen: HashSet<_> = values.iter().cloned().collect();
        values.extend(extra.into_iter().filter(|value| seen.insert(value.clone())));
    }

    fn join_annotations(value: &str) -> String {
        let mut seen = HashSet::new();
        value
            .split(';')
            .filter(|part| seen.insert(*part))
            .collect::<Vec<_>>()
            .join(";")
            .trim_start_matches(';')
            .to_string()
    }

    Ok(actual_idx)
}

/// Delete a node from the graph
pub fn delete_node(graph: &mut Graph<GeneNode, GeneEdge>, node: NodeIndex) -> Result<()> {
    if let Some(data) = graph.node_weight(node) {
        let members: Vec<_> = data.members.iter().copied().collect();
        let mut neighbours_by_member = std::collections::BTreeMap::<usize, Vec<NodeIndex>>::new();
        for edge in graph.edges(node) {
            let neighbour = edge.target();
            if neighbour != node {
                for &member in &edge.weight().members {
                    neighbours_by_member
                        .entry(member)
                        .or_default()
                        .push(neighbour);
                }
            }
        }
        for member in members {
            if let Some(neighbours) = neighbours_by_member.get_mut(&member) {
                connect_neighbours(graph, member, neighbours);
            }
        }
        graph.remove_node(node);
    }
    Ok(())
}

/// Remove a specific genome member from a node
pub fn remove_member_from_node(
    graph: &mut Graph<GeneNode, GeneEdge>,
    node: NodeIndex,
    member: usize,
) -> Result<()> {
    if graph.node_weight(node).is_none() {
        return Ok(());
    }

    if !graph[node].members.contains(&member) {
        return Ok(());
    }
    connect_member_neighbours(graph, node, member);
    let edges: Vec<_> = graph.edges(node).map(|edge| edge.id()).collect();
    for edge in edges {
        let data = graph.edge_weight_mut(edge).context("Missing graph edge")?;
        data.members.remove(&member);
        data.size = data.members.len();
        if data.members.is_empty() {
            graph.remove_edge(edge);
        }
    }

    graph[node].members.remove(&member);
    graph[node].size = graph[node].members.len();

    // Remove sequences belonging to this member
    graph[node].seq_ids.retain(|seq_id| {
        seq_id
            .split('_')
            .next()
            .and_then(|s| s.parse::<usize>().ok())
            .map(|m| m != member)
            .unwrap_or(true)
    });

    Ok(())
}

fn connect_member_neighbours(
    graph: &mut Graph<GeneNode, GeneEdge>,
    node: NodeIndex,
    member: usize,
) {
    let mut neighbours: Vec<_> = graph
        .edges(node)
        .filter(|edge| edge.weight().members.contains(&member))
        .map(|edge| edge.target())
        .filter(|&neighbour| neighbour != node)
        .collect();
    connect_neighbours(graph, member, &mut neighbours);
}

fn connect_neighbours(
    graph: &mut Graph<GeneNode, GeneEdge>,
    member: usize,
    neighbours: &mut Vec<NodeIndex>,
) {
    neighbours.sort_unstable();
    neighbours.dedup();
    for (index, &first) in neighbours.iter().enumerate() {
        for &second in &neighbours[index + 1..] {
            if let Some(edge) = graph.find_edge(first, second) {
                let data = &mut graph[edge];
                data.members.insert(member);
                data.size = data.members.len();
            } else {
                add_ordered_edge(
                    graph,
                    first,
                    second,
                    GeneEdge {
                        size: 1,
                        members: std::collections::BTreeSet::from([member]),
                        ..GeneEdge::default()
                    },
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delete_node() {
        let mut graph = Graph::default();
        let node = graph.add_node(GeneNode::default());

        assert!(graph.node_weight(node).is_some());
        delete_node(&mut graph, node).unwrap();
        assert!(graph.node_weight(node).is_none());
    }

    #[test]
    fn test_remove_member() {
        let mut graph = Graph::default();
        let mut node_data = GeneNode::default();
        node_data.members.insert(1);
        node_data.members.insert(2);
        node_data.size = 2;
        let node = graph.add_node(node_data);

        remove_member_from_node(&mut graph, node, 1).unwrap();
        assert_eq!(graph[node].members.len(), 1);
        assert!(!graph[node].members.contains(&1));
    }

    #[test]
    fn removing_member_reconnects_only_its_neighbours_and_preserves_existing_support() {
        use std::collections::BTreeSet;
        let mut graph = Graph::default();
        let center = graph.add_node(GeneNode {
            members: BTreeSet::from([0, 1]),
            size: 2,
            seq_ids: HashSet::from(["0_1".into(), "1_1".into()]),
            ..Default::default()
        });
        let left = graph.add_node(GeneNode::default());
        let right = graph.add_node(GeneNode::default());
        let other = graph.add_node(GeneNode::default());
        graph.add_edge(
            left,
            center,
            GeneEdge {
                size: 2,
                members: BTreeSet::from([0, 1]),
                ..GeneEdge::default()
            },
        );
        graph.add_edge(
            center,
            right,
            GeneEdge {
                size: 1,
                members: BTreeSet::from([0]),
                ..GeneEdge::default()
            },
        );
        graph.add_edge(
            center,
            other,
            GeneEdge {
                size: 1,
                members: BTreeSet::from([1]),
                ..GeneEdge::default()
            },
        );
        graph.add_edge(
            left,
            right,
            GeneEdge {
                size: 1,
                members: BTreeSet::from([2]),
                ..GeneEdge::default()
            },
        );
        remove_member_from_node(&mut graph, center, 0).unwrap();
        let replacement = graph.find_edge(left, right).unwrap();
        assert_eq!(graph[replacement].members, BTreeSet::from([0, 2]));
        assert_eq!(graph[replacement].size, 2);
        assert_eq!(
            graph[graph.find_edge(left, center).unwrap()].members,
            BTreeSet::from([1])
        );
        assert!(graph.find_edge(center, right).is_none());
        assert!(graph.find_edge(right, other).is_none());
        assert_eq!(graph[center].seq_ids, HashSet::from(["1_1".into()]));
    }

    #[test]
    fn deleting_node_reconnects_each_genome_separately() {
        use std::collections::BTreeSet;
        let mut graph = Graph::default();
        let center = graph.add_node(GeneNode {
            members: BTreeSet::from([0, 1]),
            size: 2,
            ..Default::default()
        });
        let a = graph.add_node(GeneNode::default());
        let b = graph.add_node(GeneNode::default());
        let c = graph.add_node(GeneNode::default());
        graph.add_edge(
            center,
            a,
            GeneEdge {
                size: 1,
                members: BTreeSet::from([0]),
                ..GeneEdge::default()
            },
        );
        graph.add_edge(
            b,
            center,
            GeneEdge {
                size: 2,
                members: BTreeSet::from([0, 1]),
                ..GeneEdge::default()
            },
        );
        graph.add_edge(
            center,
            c,
            GeneEdge {
                size: 1,
                members: BTreeSet::from([1]),
                ..GeneEdge::default()
            },
        );
        delete_node(&mut graph, center).unwrap();
        assert_eq!(
            graph[graph.find_edge(a, b).unwrap()].members,
            BTreeSet::from([0])
        );
        assert_eq!(
            graph[graph.find_edge(b, c).unwrap()].members,
            BTreeSet::from([1])
        );
        assert!(graph.find_edge(a, c).is_none());
    }

    #[test]
    fn indexed_reconnection_matches_per_member_scans() {
        use std::collections::BTreeSet;
        let mut graph = Graph::default();
        let center = graph.add_node(GeneNode {
            members: BTreeSet::from([0, 1, 2]),
            size: 3,
            ..Default::default()
        });
        let nodes: Vec<_> = (0..4)
            .map(|_| graph.add_node(GeneNode::default()))
            .collect();
        for (first, second, members) in [
            (center, nodes[0], vec![0, 1]),
            (nodes[1], center, vec![0, 2]),
            (center, nodes[1], vec![1]),
            (center, nodes[2], vec![1, 2]),
            (nodes[3], center, vec![0, 2]),
            (center, center, vec![0, 1]),
            (nodes[1], nodes[0], vec![3]),
            (nodes[0], nodes[1], vec![4]),
        ] {
            graph.add_edge(
                first,
                second,
                GeneEdge {
                    size: members.len(),
                    members: members.into_iter().collect(),
                    ..GeneEdge::default()
                },
            );
        }
        let mut expected = graph.clone();
        for member in [0, 1, 2] {
            connect_member_neighbours(&mut expected, center, member);
        }
        expected.remove_node(center);
        delete_node(&mut graph, center).unwrap();
        let snapshot = |graph: &Graph<GeneNode, GeneEdge>| {
            graph
                .edge_indices()
                .map(|edge| {
                    (
                        edge,
                        graph.edge_endpoints(edge),
                        graph[edge].size,
                        graph[edge].members.clone(),
                    )
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(snapshot(&graph), snapshot(&expected));
    }
}
