//! CD-HIT wrapper module
//!
//! Provides functions to invoke cd-hit and cd-hit-est binaries

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::process::Stdio;

pub(crate) struct CdhitConfig {
    pub identity: f64,
    pub length_ratio: f64,
    pub threads: usize,
    pub dna: bool,
    pub accurate: bool,
    pub quiet: bool,
}

pub(crate) fn run_cluster(input: &str, output: &str, config: &CdhitConfig) -> Result<()> {
    let binary = if config.dna { "cd-hit-est" } else { "cd-hit" };
    let mut command = Command::new(binary);
    command
        .args(["-i", input, "-o", output, "-c"])
        .arg(config.identity.to_string())
        .arg("-s")
        .arg(config.length_ratio.to_string())
        .args(["-aL", "0", "-AL", "99999999", "-aS", "0", "-AS", "99999999"])
        .args(["-M", "0", "-d", "999", "-T"])
        .arg(config.threads.to_string());
    if config.accurate {
        command.args(["-g", "1", "-n"]);
        command.arg(if config.dna {
            "6"
        } else {
            cdhit_word_length(config.identity)
        });
    } else if config.dna {
        command.args(["-n", "7"]);
    }
    if config.dna {
        command.args(["-r", "1", "-mask", "NX"]);
    }
    if config.quiet {
        command.stdout(Stdio::null());
    }
    let status = command
        .status()
        .with_context(|| format!("Failed to run {}", binary))?;
    if !status.success() {
        anyhow::bail!(
            "{} failed at identity {} with status {}",
            binary,
            config.identity,
            status
        );
    }
    Ok(())
}

/// Check CD-HIT version
pub fn check_cdhit_version(binary: &str) -> Result<String> {
    let output = Command::new(binary)
        .arg("-h")
        .output()
        .with_context(|| format!("Failed to run {}. Is it installed?", binary))?;

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    Ok(format!("{}{}", stdout, stderr))
}

/// Run CD-HIT on protein sequences
pub fn run_cdhit(
    input_file: &str,
    output_file: &str,
    id: f64,
    s: f64,
    n_cpu: usize,
    _outdir: &str,
) -> Result<()> {
    run_cluster(
        input_file,
        output_file,
        &CdhitConfig {
            identity: id,
            length_ratio: s,
            threads: n_cpu,
            dna: false,
            accurate: true,
            quiet: true,
        },
    )
}

fn cdhit_word_length(identity: f64) -> &'static str {
    if identity >= 0.7 {
        "5"
    } else if identity >= 0.6 {
        "3"
    } else {
        "2"
    }
}

/// Run CD-HIT-EST on nucleotide sequences
pub fn run_cdhit_est(
    input_file: &str,
    output_file: &str,
    id: f64,
    s: f64,
    strand: i32,
    n_cpu: usize,
) -> Result<()> {
    let status = Command::new("cd-hit-est")
        .arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .arg("-c")
        .arg(format!("{:.2}", id))
        .arg("-s")
        .arg(format!("{:.2}", s))
        .arg("-r")
        .arg(strand.to_string())
        .arg("-T")
        .arg(n_cpu.to_string())
        .arg("-M")
        .arg("0")
        .arg("-g")
        .arg("1")
        .arg("-d")
        .arg("0")
        .status()
        .context("Failed to run cd-hit-est")?;

    if !status.success() {
        anyhow::bail!("cd-hit-est failed with exit code: {:?}", status.code());
    }

    Ok(())
}

/// Parse CD-HIT .clstr file
pub fn parse_cdhit_clusters(clstr_file: &str) -> Result<HashMap<String, Vec<String>>> {
    let file = File::open(clstr_file)?;
    let reader = BufReader::new(file);

    let mut clusters = HashMap::new();
    let mut current_cluster = String::new();
    let mut current_members = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_cluster.is_empty() {
                clusters.insert(current_cluster.clone(), current_members.clone());
            }
            current_cluster = line[1..].trim().to_string();
            current_members = Vec::new();
        } else if line.contains('>') {
            let parts: Vec<&str> = line.split('>').collect();
            if parts.len() >= 2 {
                let seq_id = parts[1].split("...").next().unwrap_or("").trim();
                current_members.push(seq_id.to_string());
            }
        }
    }

    if !current_cluster.is_empty() {
        clusters.insert(current_cluster, current_members);
    }

    Ok(clusters)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cdhit_word_lengths_support_low_thresholds() {
        assert_eq!(cdhit_word_length(0.99), "5");
        assert_eq!(cdhit_word_length(0.70), "5");
        assert_eq!(cdhit_word_length(0.60), "3");
        assert_eq!(cdhit_word_length(0.50), "2");
    }

    #[test]
    fn test_version_check() {
        // This will fail if cd-hit is not installed, which is expected
        let _ = check_cdhit_version("cd-hit");
    }
}
