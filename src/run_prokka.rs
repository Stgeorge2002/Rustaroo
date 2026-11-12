//! Run Prokka gene annotation
//!

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;
use which::which;

/// Run Prokka with modified prodigal settings
pub fn run_prokka_mod(
    input_file: &str,
    out_folder: &str,
    train_file: &str,
    force: bool,
    add_cmds: Option<&str>,
) -> Result<()> {
    let input_path = Path::new(input_file);
    let prefix = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow::anyhow!("Invalid input file name"))?;

    let path_to_prodigal = which("prodigal").context("Prodigal not found in PATH")?;

    // Create temporary directory for modified prodigal wrapper
    let temp_dir = TempDir::new_in(out_folder)?;
    let temp_path = temp_dir.path();

    // Create wrapper script for prodigal that includes training file
    let wrapper_script = temp_path.join("prodigal");
    let mut script = fs::File::create(&wrapper_script)?;
    writeln!(script, "#!/bin/bash")?;
    writeln!(script, "(>&2 echo 'running prokka mod!')")?;
    writeln!(
        script,
        "{} -t {} \"$@\"",
        path_to_prodigal.display(),
        Path::new(train_file).canonicalize()?.display()
    )?;
    script.sync_all()?;

    // Make wrapper executable
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = fs::metadata(&wrapper_script)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&wrapper_script, perms)?;
    }

    // Build prokka command
    let mut cmd = Command::new("prokka");

    if let Some(additional) = add_cmds {
        // Parse additional commands carefully
        for arg in additional.split_whitespace() {
            cmd.arg(arg);
        }
    }

    if force {
        cmd.arg("--force");
    }

    cmd.arg("--cpus").arg("1");
    cmd.arg("--outdir").arg(format!("{}{}", out_folder, prefix));
    cmd.arg("--prefix").arg(prefix);
    cmd.arg(input_file);

    // Set PATH to include our temporary directory first
    let path_var = std::env::var("PATH").unwrap_or_default();
    let new_path = format!("{}:{}", temp_path.display(), path_var);
    cmd.env("PATH", new_path);
    cmd.env("TMPDIR", temp_path);

    // Redirect output to log file
    let log_file = format!("{}{}_prokka.log", out_folder, prefix);
    let log = fs::File::create(&log_file)?;
    cmd.stdout(log.try_clone()?);
    cmd.stderr(log);

    // Run prokka
    let status = cmd.status().context("Failed to run prokka")?;

    if !status.success() {
        anyhow::bail!("Prokka failed with status: {}", status);
    }

    // Check prokka completed successfully
    let log_contents = fs::read_to_string(&log_file)?;
    let lines: Vec<&str> = log_contents.lines().collect();
    if lines.len() < 6 {
        anyhow::bail!("Prokka log file too short, likely failed");
    }
    let last_lines = &lines[lines.len().saturating_sub(6)..];
    if !last_lines
        .iter()
        .any(|line| line.contains("Annotation finished successfully."))
    {
        anyhow::bail!("Prokka did not execute successfully!");
    }

    Ok(())
}

/// Run Prokka on multiple genomes in parallel
pub fn run_prokka_parallel(
    input_files: &[String],
    output_dir: &str,
    training_file: &str,
    force: bool,
    add_cmds: Option<&str>,
    n_cpu: usize,
) -> Result<()> {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_cpu).build()?;

    pool.install(|| {
        input_files.par_iter().try_for_each(|input_file| {
            run_prokka_mod(input_file, output_dir, training_file, force, add_cmds)
        })
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_prokka_mod_setup() {
        // Test would require prokka installed
        // Just verify functions exist
    }
}
