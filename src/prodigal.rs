//! Prodigal gene prediction wrapper
//!

use anyhow::{Context, Result};
use rand::seq::SliceRandom;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::Command;
use tempfile::NamedTempFile;

/// Run Prodigal on the command line
#[allow(clippy::too_many_arguments)]
pub fn run_prodigal(
    trans_file: &str,
    nuc_file: &str,
    input_file: &str,
    output_file: &str,
    mode: &str,
    shine_dalgarno: bool,
    mask_sequence: bool,
    tr_table: u8,
    closed_ends: bool,
    output_type: &str,
    quiet: bool,
    start_file: Option<&str>,
    training_file: Option<&str>,
) -> Result<()> {
    let mut cmd = Command::new("prodigal");

    cmd.arg("-a")
        .arg(trans_file)
        .arg("-d")
        .arg(nuc_file)
        .arg("-f")
        .arg(output_type)
        .arg("-g")
        .arg(tr_table.to_string())
        .arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .arg("-p")
        .arg(mode);

    if closed_ends {
        cmd.arg("-c");
    }

    if mask_sequence {
        cmd.arg("-m");
    }

    if !shine_dalgarno {
        cmd.arg("-n");
    }

    if let Some(start) = start_file {
        cmd.arg("-s").arg(start);
    }

    if let Some(training) = training_file {
        cmd.arg("-t").arg(training);
    }

    if quiet {
        cmd.arg("-q");
    } else {
        println!("running cmd: prodigal {:?}", cmd);
    }

    let status = cmd.status().context("Failed to run prodigal")?;

    if !status.success() {
        anyhow::bail!("Prodigal failed with status: {}", status);
    }

    Ok(())
}

/// Train Prodigal using a sample of input genomes
pub fn train_prodigal(
    input_files: &[String],
    n_samples: usize,
    force: bool,
    outdir: &str,
) -> Result<()> {
    let training_file = format!("{}prodigal_training.txt", outdir);

    if Path::new(&training_file).exists() {
        if force {
            fs::remove_file(&training_file)?;
        } else {
            // Training file already exists
            return Ok(());
        }
    }

    // Randomly sample training genomes
    let samples: Vec<_> = if input_files.len() < n_samples {
        input_files.to_vec()
    } else {
        let mut rng = rand::thread_rng();
        let mut files = input_files.to_vec();
        files.shuffle(&mut rng);
        files.into_iter().take(n_samples).collect()
    };

    // Combine training files into one fasta
    let mut temp_file = NamedTempFile::new_in(outdir)?;
    for infile in &samples {
        let contents = fs::read_to_string(infile)?;
        writeln!(temp_file, "{}", contents.trim())?;
    }
    temp_file.flush()?;

    // Run prodigal
    run_prodigal(
        "/dev/null",
        "/dev/null",
        temp_file.path().to_str().unwrap(),
        "/dev/null",
        "single",
        true,
        false,
        11,
        false,
        "gff",
        false,
        None,
        Some(&training_file),
    )?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prodigal_command_construction() {
        // This test would require prodigal to be installed
        // Just verify the function exists and can be called
        let _result = run_prodigal(
            "test.faa",
            "test.fna",
            "input.fasta",
            "output.gff",
            "single",
            true,
            false,
            11,
            false,
            "gff",
            true,
            None,
            None,
        );
        // Will fail without prodigal installed, but that's expected
    }
}
