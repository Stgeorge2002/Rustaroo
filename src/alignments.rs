//! Generate multiple sequence alignments for core genes
//!

use crate::network::Graph;
use crate::network::{GeneEdge, GeneNode};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::Write;
use std::process::Command;

/// Generate alignments using MAFFT or PRANK
pub fn generate_core_gene_alignments(
    graph: &Graph<GeneNode, GeneEdge>,
    output_dir: &str,
    aligner: &str,
    n_cpu: usize,
    core_threshold: f64,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        println!("Generating core gene alignments with {}...", aligner);
    }

    let n_isolates = graph
        .node_indices()
        .map(|n| graph[n].members.len())
        .max()
        .unwrap_or(0);

    let core_cutoff = (n_isolates as f64 * core_threshold) as usize;

    // Find core genes
    let core_genes: Vec<_> = graph
        .node_indices()
        .filter(|&n| graph[n].size >= core_cutoff)
        .collect();

    if !quiet {
        println!("Found {} core genes", core_genes.len());
    }

    // Generate alignment for each core gene
    for (i, &node) in core_genes.iter().enumerate() {
        let sequences = &graph[node].dna;
        if sequences.len() < 2 {
            continue;
        }

        let input_file = format!("{}/core_gene_{}.fasta", output_dir, i);
        let output_file = format!("{}/core_gene_{}.aln", output_dir, i);

        // Write sequences to temp file
        let mut file = File::create(&input_file)?;
        for (j, seq) in sequences.iter().enumerate() {
            writeln!(file, ">seq_{}", j)?;
            writeln!(file, "{}", seq)?;
        }
        drop(file);

        // Run aligner
        match aligner {
            "mafft" => run_mafft(&input_file, &output_file, n_cpu)?,
            "prank" => run_prank(&input_file, &output_file)?,
            _ => anyhow::bail!("Unknown aligner: {}", aligner),
        }
    }

    if !quiet {
        println!("Alignment generation complete");
    }

    Ok(())
}

/// Run MAFFT
fn run_mafft(input: &str, output: &str, n_cpu: usize) -> Result<()> {
    let status = Command::new("mafft")
        .arg("--auto")
        .arg("--thread")
        .arg(n_cpu.to_string())
        .arg(input)
        .stdout(File::create(output)?)
        .status()
        .context("Failed to run mafft")?;

    if !status.success() {
        anyhow::bail!("mafft failed");
    }

    Ok(())
}

/// Run PRANK
fn run_prank(input: &str, output: &str) -> Result<()> {
    let status = Command::new("prank")
        .arg("-d=")
        .arg(input)
        .arg(format!("-o={}", output))
        .status()
        .context("Failed to run prank")?;

    if !status.success() {
        anyhow::bail!("prank failed");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_alignment_module() {
        // Test requires actual aligner binaries
    }
}
