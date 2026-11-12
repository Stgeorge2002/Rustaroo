//! Estimate parameters for genome presence inference
//!
//! Port of scripts/estimate_img_params.py

use anyhow::Result;

/// Estimate parameters for inferring missing genomes
pub fn estimate_img_parameters(
    gene_presence_absence: &str,
    output_file: &str,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        println!("Estimating IMG parameters...");
    }

    use std::fs::File;
    use std::io::{BufRead, BufReader, Write};

    let file = File::open(gene_presence_absence)?;
    let reader = BufReader::new(file);

    let mut n_genes = 0;
    let mut n_samples = 0;
    let mut presence_counts = Vec::new();

    for (i, line_result) in reader.lines().enumerate() {
        let line = line_result?;

        if i == 0 {
            let fields: Vec<&str> = line.split('\t').collect();
            n_samples = fields.len().saturating_sub(1);
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() > 1 {
            n_genes += 1;
            let present = fields[1..].iter().filter(|&&f| f != "0").count();
            presence_counts.push(present);
        }
    }

    let mean_presence = if !presence_counts.is_empty() {
        presence_counts.iter().sum::<usize>() as f64 / presence_counts.len() as f64
    } else {
        0.0
    };

    let mut output = File::create(output_file)?;
    writeln!(output, "n_genes: {}", n_genes)?;
    writeln!(output, "n_samples: {}", n_samples)?;
    writeln!(output, "mean_presence: {:.2}", mean_presence)?;

    if !quiet {
        println!("Estimated parameters written to {}", output_file);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_estimate_img() {
        // Test requires actual PA file
    }
}
