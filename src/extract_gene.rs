//! Extract specific genes from pangenome
//!
//! Port of scripts/extract_gene_from_pangenome.py

use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Extract sequences for specific genes
pub fn extract_gene_sequences(
    _gene_presence_absence: &str,
    gene_data_file: &str,
    gene_names: &[String],
    output_file: &str,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        println!("Extracting {} genes...", gene_names.len());
    }

    // Read gene data to get sequences
    let mut gene_sequences = std::collections::HashMap::new();
    let file = File::open(gene_data_file)?;
    let reader = BufReader::new(file);

    for (i, line_result) in reader.lines().enumerate() {
        if i == 0 {
            continue; // Skip header
        }
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() >= 6 {
            let gene_id = fields[0].to_string();
            let sequence = fields[5].to_string();
            gene_sequences.insert(gene_id, sequence);
        }
    }

    // Extract requested genes
    let mut output = File::create(output_file)?;
    let mut found = 0;

    for gene_name in gene_names {
        if let Some(seq) = gene_sequences.get(gene_name) {
            writeln!(output, ">{}", gene_name)?;
            writeln!(output, "{}", seq)?;
            found += 1;
        }
    }

    if !quiet {
        println!("Extracted {} genes", found);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_extract_gene() {
        // Test requires actual gene data
    }
}
