//! Integrate new genomes into existing pangenome
//!

use crate::network::Graph;
use crate::network::{GeneEdge, GeneNode};
use anyhow::Result;
use std::path::Path;

/// Integrate new isolates into an existing pangenome graph
pub fn integrate_new_genes(
    _graph: &mut Graph<GeneNode, GeneEdge>,
    new_gff_files: &[String],
    _output_dir: &str,
    _n_cpu: usize,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        println!("Integrating {} new isolates...", new_gff_files.len());
    }

    // Process new GFF files
    for gff_file in new_gff_files {
        if !Path::new(gff_file).exists() {
            anyhow::bail!("GFF file not found: {}", gff_file);
        }

        // TODO: Parse GFF, extract genes, cluster against existing nodes, add to graph
        if !quiet {
            println!("Processing {}...", gff_file);
        }
    }

    if !quiet {
        println!("Integration complete");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integrate() {
        let mut graph: Graph<GeneNode, GeneEdge> = Graph::default();
        let result = integrate_new_genes(&mut graph, &[], "", 1, true);
        assert!(result.is_ok());
    }
}
