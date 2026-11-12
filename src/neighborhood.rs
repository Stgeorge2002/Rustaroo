//! Analyze gene neighborhoods
//!
//! Port of scripts/gene_neighbourhood.py

use crate::network::Graph;
use crate::network::{GeneEdge, GeneNode};
use anyhow::Result;

/// Extract gene neighborhood information
pub fn extract_gene_neighborhood(
    graph: &Graph<GeneNode, GeneEdge>,
    gene_name: &str,
    radius: usize,
    output_file: &str,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        println!("Extracting neighborhood for gene: {}", gene_name);
    }

    // Find the gene node
    let mut target_node = None;
    for node_idx in graph.node_indices() {
        if graph[node_idx].annotation.contains(gene_name) {
            target_node = Some(node_idx);
            break;
        }
    }

    let target = match target_node {
        Some(n) => n,
        None => anyhow::bail!("Gene {} not found", gene_name),
    };

    // BFS to find neighbors within radius
    let mut neighbors = Vec::new();
    let mut visited = std::collections::HashSet::new();
    let mut queue = std::collections::VecDeque::new();

    queue.push_back((target, 0));
    visited.insert(target);

    while let Some((node, dist)) = queue.pop_front() {
        if dist <= radius {
            neighbors.push((node, dist));

            if dist < radius {
                for neighbor in graph.neighbors(node) {
                    if !visited.contains(&neighbor) {
                        visited.insert(neighbor);
                        queue.push_back((neighbor, dist + 1));
                    }
                }
            }
        }
    }

    if !quiet {
        println!(
            "Found {} neighbors within radius {}",
            neighbors.len(),
            radius
        );
    }

    // Write neighborhood to file
    use std::fs::File;
    use std::io::Write;
    let mut file = File::create(output_file)?;
    writeln!(file, "gene,distance,annotation")?;
    for (node, dist) in neighbors {
        writeln!(
            file,
            "{},{},{}",
            graph[node].annotation, dist, graph[node].description
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_neighborhood() {
        // Test requires actual graph
    }
}
