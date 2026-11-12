//! Node merging utilities for graph manipulation
//!

use crate::network::Graph;
use crate::network::{GeneEdge, GeneNode};
use anyhow::{Context, Result};
use petgraph::graph::NodeIndex;
use std::collections::HashSet;

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
    let mut neighbor_edges: std::collections::HashMap<NodeIndex, GeneEdge> =
        std::collections::HashMap::new();

    for &node_idx in &cluster {
        if graph.node_weight(node_idx).is_none() {
            continue;
        }

        use petgraph::visit::EdgeRef;
        for edge in graph.edges(node_idx) {
            let neighbor = if edge.source() == node_idx {
                edge.target()
            } else {
                edge.source()
            };

            // Skip edges within the cluster
            if cluster_set.contains(&neighbor) {
                continue;
            }

            neighbor_edges
                .entry(neighbor)
                .and_modify(|e| {
                    e.members.extend(edge.weight().members.iter());
                    e.size = e.members.len();
                })
                .or_insert_with(|| edge.weight().clone());
        }
    }

    let actual_idx = graph.add_node(merged);

    // Add edges to neighbors
    for (neighbor, edge_data) in neighbor_edges {
        if graph.node_weight(neighbor).is_some() {
            graph.add_edge(actual_idx, neighbor, edge_data);
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
    if graph.node_weight(node).is_some() {
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

    graph[node].members.remove(&member);
    if graph[node].size > 0 {
        graph[node].size -= 1;
    }

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
}
