use anyhow::{ensure, Result};
use std::fs;
use std::io;
use std::path::Path;
use std::process::Command;

fn rustaroo(out: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_rustaroo"));
    command
        .arg("--out-dir")
        .arg(out)
        .args(["--clean-mode", "strict", "--skip-core-aln"]);
    command
}

#[test]
#[ignore = "requires cd-hit and cd-hit-est on PATH"]
fn minimal_retains_only_the_normal_reference_matrix_and_log() -> Result<()> {
    let temp = tempfile::tempdir()?;
    let dna = format!("ATG{}TAA", "GCTGATTTTCAACGT".repeat(12));
    let mut inputs = Vec::new();
    for isolate in 0..3 {
        let input = temp.path().join(format!("isolate_{isolate}.gff"));
        fs::write(
            &input,
            format!(
                "##gff-version 3\ncontig\tProkka\tCDS\t1\t{}\t.\t+\t0\tID=gene_{isolate};gene=example;product=Example\n##FASTA\n>contig\n{dna}\n",
                dna.len()
            ),
        )?;
        inputs.push(input);
    }

    let full = temp.path().join("full");
    let minimal = temp.path().join("minimal");
    for (out, flags) in [(&full, &[][..]), (&minimal, &["--minimal", "--quiet"][..])] {
        let output = rustaroo(out).args(flags).arg("-i").args(&inputs).output()?;
        ensure!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        if out == &minimal {
            assert!(output.stderr.is_empty());
        }
    }

    let mut files = fs::read_dir(&minimal)?
        .map(|entry| entry.map(|entry| entry.file_name()))
        .collect::<io::Result<Vec<_>>>()?;
    files.sort();
    assert_eq!(
        files,
        [
            "gene_presence_absence.csv",
            "pan_genome_reference.fa",
            "rustaroo.log"
        ]
    );
    for file in ["gene_presence_absence.csv", "pan_genome_reference.fa"] {
        let expected = fs::read(full.join(file))?;
        assert!(!expected.is_empty(), "{file}");
        assert_eq!(fs::read(minimal.join(file))?, expected, "{file}");
    }
    let log = fs::read_to_string(minimal.join("rustaroo.log"))?;
    assert!(log.contains("Pipeline total:"));
    assert!(log.contains("Rustaroo completed successfully!"));
    Ok(())
}

#[test]
fn minimal_keeps_a_log_when_the_run_fails() -> Result<()> {
    let temp = tempfile::tempdir()?;
    let out = temp.path().join("failed");
    let output = rustaroo(&out)
        .args(["--minimal", "--quiet", "-i"])
        .arg(temp.path().join("missing.gff"))
        .output()?;
    assert!(!output.status.success());
    assert!(fs::read_to_string(out.join("rustaroo.log"))?.contains("Input file does not exist"));
    assert_eq!(fs::read_dir(out)?.count(), 1);
    Ok(())
}
