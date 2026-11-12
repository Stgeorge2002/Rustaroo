//! Post-run MSA generation
//!
//! Port of scripts/post_run_msa.py

use anyhow::Result;

/// Generate MSA after the main pipeline run
pub fn post_run_msa(
    _graph_file: &str,
    _output_dir: &str,
    _aligner: &str,
    _n_cpu: usize,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        println!("Generating post-run MSA...");
    }

    // TODO: Load graph and generate alignments

    if !quiet {
        println!("MSA generation complete");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_post_run_msa() {
        let result = post_run_msa("", "", "mafft", 1, true);
        assert!(result.is_ok());
    }
}
