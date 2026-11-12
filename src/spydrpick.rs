//! SpydrPick input file generation
//!
//! Port of scripts/generate_spydrpick_input.py

use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Generate SpydrPick input from presence/absence matrix
pub fn generate_spydrpick_input(
    gene_presence_absence: &str,
    output_file: &str,
    min_variant_support: usize,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        println!("Generating SpydrPick input...");
    }

    let file = File::open(gene_presence_absence)?;
    let reader = BufReader::new(file);
    let mut output = File::create(output_file)?;

    let mut gene_count = 0;

    for (i, line_result) in reader.lines().enumerate() {
        let line = line_result?;

        if i == 0 {
            // Header - write isolate names
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() > 1 {
                writeln!(output, "{}", fields[1..].join("\t"))?;
            }
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }

        // Count variants
        let present = fields[1..].iter().filter(|&&f| f != "0").count();
        let absent = fields[1..].iter().filter(|&&f| f == "0").count();

        // Only include genes with sufficient variation
        if present >= min_variant_support && absent >= min_variant_support {
            writeln!(output, "{}", fields[1..].join("\t"))?;
            gene_count += 1;
        }
    }

    if !quiet {
        println!("Generated SpydrPick input with {} genes", gene_count);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_spydrpick() {
        // Test requires actual PA matrix
    }
}
