//! Filter presence/absence matrix
//!
//! Port of scripts/filter_presence_absence.py

use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Filter presence/absence matrix by various criteria
pub fn filter_presence_absence(
    input_file: &str,
    output_file: &str,
    min_samples: Option<usize>,
    max_samples: Option<usize>,
    quiet: bool,
) -> Result<()> {
    let file = File::open(input_file)?;
    let reader = BufReader::new(file);
    let mut output = File::create(output_file)?;

    let mut kept = 0;
    let mut filtered = 0;

    for (i, line_result) in reader.lines().enumerate() {
        let line = line_result?;

        if i == 0 {
            // Header
            writeln!(output, "{}", line)?;
            continue;
        }

        // Count presence
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }

        let present_count = fields[1..]
            .iter()
            .filter(|&&f| f != "0" && !f.is_empty())
            .count();

        let mut keep = true;
        if let Some(min) = min_samples {
            if present_count < min {
                keep = false;
            }
        }
        if let Some(max) = max_samples {
            if present_count > max {
                keep = false;
            }
        }

        if keep {
            writeln!(output, "{}", line)?;
            kept += 1;
        } else {
            filtered += 1;
        }
    }

    if !quiet {
        println!("Kept {} genes, filtered {} genes", kept, filtered);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_filter_pa() {
        // Test requires actual PA file
    }
}
