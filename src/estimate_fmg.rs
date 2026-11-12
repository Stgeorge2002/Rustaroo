//! Estimate parameters for finding many genes
//!
//! Port of scripts/estimate_many_genes_params.py

use anyhow::Result;

/// Estimate parameters for genes present in many genomes
pub fn estimate_many_genes_parameters(
    gene_presence_absence: &str,
    output_file: &str,
    threshold: f64,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        println!("Estimating many-genes parameters...");
    }

    use std::fs::File;
    use std::io::{BufRead, BufReader, Write};

    let file = File::open(gene_presence_absence)?;
    let reader = BufReader::new(file);

    let mut n_samples = 0;
    let mut high_frequency_genes = Vec::new();

    for (i, line_result) in reader.lines().enumerate() {
        let line = line_result?;

        if i == 0 {
            let fields: Vec<&str> = line.split('\t').collect();
            n_samples = fields.len().saturating_sub(1);
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() > 1 {
            let gene_name = fields[0];
            let present = fields[1..].iter().filter(|&&f| f != "0").count();
            let freq = present as f64 / n_samples as f64;

            if freq >= threshold {
                high_frequency_genes.push((gene_name.to_string(), freq));
            }
        }
    }

    let mut output = File::create(output_file)?;
    writeln!(output, "threshold: {}", threshold)?;
    writeln!(
        output,
        "high_frequency_genes: {}",
        high_frequency_genes.len()
    )?;
    writeln!(output, "n_samples: {}", n_samples)?;

    if !quiet {
        println!("Found {} high-frequency genes", high_frequency_genes.len());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_estimate_fmg() {
        // Test requires actual PA file
    }
}
