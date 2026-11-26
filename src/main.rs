//! Main Rustaroo CLI
//!

use anyhow::{Context, Result};
use clap::Parser;
use rayon::prelude::*;
use rustaroo::*;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq)]
enum CleanMode {
    Strict,
    Moderate,
    Sensitive,
}

impl std::str::FromStr for CleanMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "strict" => Ok(Self::Strict),
            "moderate" => Ok(Self::Moderate),
            "sensitive" => Ok(Self::Sensitive),
            _ => Err(format!("Unknown mode: {}", s)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum RefindMode {
    Default,
    Strict,
    Off,
}

impl std::str::FromStr for RefindMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "default" => Ok(Self::Default),
            "strict" => Ok(Self::Strict),
            "off" => Ok(Self::Off),
            _ => Err(format!("Unknown refind mode: {}", s)),
        }
    }
}

#[derive(Parser, Debug)]
#[command(name = "rustaroo")]
#[command(about = "An updated pipeline for pangenome investigation", long_about = None)]
#[command(version = rustaroo::VERSION)]
struct Args {
    /// Input GFF3 files (usually output from Prokka)
    #[arg(short = 'i', long, required = true, num_args = 1..)]
    input: Vec<PathBuf>,

    /// Output directory
    #[arg(short = 'o', long)]
    out_dir: PathBuf,

    /// Keep only rustaroo.log, pan_genome_reference.fa, and gene_presence_absence.csv
    #[arg(long, conflicts_with = "keep_temp")]
    minimal: bool,

    /// Stringency mode (strict/moderate/sensitive)
    #[arg(long = "clean-mode")]
    mode: CleanMode,

    /// Remove annotations that don't conform to Prokka format
    #[arg(long)]
    remove_invalid_genes: bool,

    /// Sequence identity threshold
    #[arg(short = 'c', long)]
    threshold: Option<f64>,

    /// Protein family sequence identity threshold
    #[arg(short = 'f', long)]
    family_threshold: Option<f64>,

    /// Length difference cutoff
    #[arg(long)]
    len_dif_percent: Option<f64>,

    /// Length difference cutoff at gene family level
    #[arg(long, default_value = "0.0")]
    family_len_dif_percent: f64,

    /// Don't split paralogs
    #[arg(long)]
    merge_paralogs: bool,

    /// Search radius for refinding genes (nucleotides)
    #[arg(long, default_value = "5000")]
    search_radius: usize,

    /// Proportion of gene required to match when refinding
    #[arg(long, default_value = "0.2")]
    refind_prop_match: f64,

    /// Refind mode
    #[arg(long = "refind-mode", default_value = "default")]
    refind_mode: RefindMode,

    /// Break equal-support refinding overlaps by node creation order, newest first
    #[arg(long)]
    refind_creation_order: bool,

    /// Minimum cluster size to keep trailing genes
    #[arg(long)]
    min_trailing_support: Option<usize>,

    /// Number of times to trim low support trailing nodes
    #[arg(long)]
    trailing_recursive: Option<usize>,

    /// Minimum edge support threshold
    #[arg(long)]
    edge_support_threshold: Option<usize>,

    /// Length outlier support proportion
    #[arg(long, default_value = "0.1")]
    length_outlier_support_proportion: f64,

    /// Remove genes by consensus
    #[arg(long)]
    remove_by_consensus: Option<bool>,

    /// Number of CPU threads to use
    #[arg(short = 't', long, default_value = "1")]
    threads: usize,

    /// Verbose output
    #[arg(long)]
    verbose: bool,

    /// Quiet mode
    #[arg(long)]
    quiet: bool,

    /// Aligner for MSA (mafft/prank)
    #[arg(short = 'a', long, default_value = "mafft")]
    aligner: String,

    /// Skip core gene alignment
    #[arg(long)]
    skip_core_aln: bool,

    /// Genetic code table (default: 11)
    #[arg(long, default_value = "11")]
    codon_table: usize,

    /// Keep temporary files
    #[arg(long)]
    keep_temp: bool,
}

struct RunLog {
    file: std::fs::File,
    console: bool,
    error: Arc<parking_lot::Mutex<Option<io::Error>>>,
}

impl Write for RunLog {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if let Err(error) = self.file.write_all(buf) {
            // env_logger cannot return sink errors to the caller.
            self.error.lock().get_or_insert(error);
        }
        if self.console {
            io::stderr().write_all(buf)?;
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.threads == 0 {
        anyhow::bail!("--threads must be greater than zero");
    }

    let log_path = args.out_dir.join("rustaroo.log");
    let log_error = Arc::new(parking_lot::Mutex::new(None));
    if args.minimal || args.verbose || !args.quiet {
        let mut logger = env_logger::Builder::from_default_env();
        logger.filter_level(if args.verbose {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Info
        });
        if args.minimal {
            std::fs::create_dir_all(&args.out_dir)?;
            logger
                .write_style(env_logger::WriteStyle::Never)
                .target(env_logger::Target::Pipe(Box::new(RunLog {
                    file: std::fs::File::create(&log_path)
                        .with_context(|| format!("Creating {}", log_path.display()))?,
                    console: args.verbose || !args.quiet,
                    error: Arc::clone(&log_error),
                })));
        }
        logger.init();
    }

    let result = run(&args);
    if args.minimal {
        if let Err(error) = &result {
            log::error!("{error:#}");
        }
        if let Some(error) = log_error.lock().take() {
            if let Err(pipeline_error) = &result {
                return Err(error).with_context(|| {
                    format!(
                        "Writing {} after run failed: {pipeline_error:#}",
                        log_path.display()
                    )
                });
            }
            return Err(error).with_context(|| format!("Writing {}", log_path.display()));
        }
    }
    result
}

fn run(args: &Args) -> Result<()> {
    log::info!("Rustaroo v{}", rustaroo::VERSION);
    log::info!("Mode: {:?}", args.mode);

    // Validate input files
    for input_file in &args.input {
        if !input_file.exists() {
            anyhow::bail!("Input file does not exist: {}", input_file.display());
        }
    }

    std::fs::create_dir_all(&args.out_dir)?;
    let work_dir = if args.minimal {
        Some(
            tempfile::Builder::new()
                .prefix(".rustaroo-")
                .tempdir_in(&args.out_dir)?,
        )
    } else {
        None
    };
    let work_path = work_dir
        .as_ref()
        .map_or(args.out_dir.as_path(), |dir| dir.path());

    // Convert mode to string for compatibility
    let mode_str = match args.mode {
        CleanMode::Strict => "strict",
        CleanMode::Moderate => "moderate",
        CleanMode::Sensitive => "sensitive",
    };

    // Set up rustaroo arguments
    let input_files: Vec<String> = args
        .input
        .iter()
        .map(|p| p.to_string_lossy().to_string())
        .collect();

    let mut rustaroo_args = args::RustarooArgs::new(mode_str.to_string(), input_files);
    rustaroo_args.id = args.threshold;
    rustaroo_args.family_threshold = args.family_threshold;
    rustaroo_args.len_dif_percent = args.len_dif_percent;
    rustaroo_args.min_trailing_support = args.min_trailing_support;
    rustaroo_args.trailing_recursive = args.trailing_recursive;
    rustaroo_args.edge_support_threshold = args.edge_support_threshold;
    rustaroo_args.remove_by_consensus = args.remove_by_consensus;

    // Set default args based on mode
    let rustaroo_args = args::set_default_args(rustaroo_args);

    log::info!("Processing {} input files", rustaroo_args.input_files.len());

    // Main pipeline execution
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    pool.install(|| run_rustaroo_pipeline(&rustaroo_args, args, work_path))?;
    if let Some(work_dir) = work_dir {
        for name in ["pan_genome_reference.fa", "gene_presence_absence.csv"] {
            let destination = args.out_dir.join(name);
            std::fs::rename(work_dir.path().join(name), &destination)
                .with_context(|| format!("Publishing {}", destination.display()))?;
        }
        work_dir
            .close()
            .context("Removing minimal-mode working files")?;
    }

    log::info!("Rustaroo completed successfully!");
    log::info!("Output written to: {}", args.out_dir.display());

    Ok(())
}

fn run_rustaroo_pipeline(
    rustaroo_args: &args::RustarooArgs,
    args: &Args,
    outdir_path: &Path,
) -> Result<()> {
    let pipeline_started = Instant::now();
    let mut timings = Vec::new();
    let outdir = outdir_path.to_string_lossy();

    // Step 1: Process input GFF files
    log::info!("Step 1: Processing input GFF files...");
    let stage_started = Instant::now();
    let all_genes: Vec<usize> = rustaroo_args
        .input_files
        .par_iter()
        .enumerate()
        .map(|(i, input_file)| {
            let proteins =
                prokka::get_gene_list(input_file, i, args.remove_invalid_genes, args.codon_table)?;

            // Write intermediate output files
            let prot_file = outdir_path.join(format!("input_{}.prot.fa", i));
            let dna_file = outdir_path.join(format!("input_{}.dna.fa", i));
            let csv_file = outdir_path.join(format!("input_{}.csv", i));

            if !proteins.is_empty() {
                prokka::output_files(
                    &std::collections::HashMap::new(),
                    &proteins,
                    prot_file.to_str().unwrap(),
                    dna_file.to_str().unwrap(),
                    csv_file.to_str().unwrap(),
                    input_file,
                )?;
            } else {
                anyhow::bail!("No valid coding sequences in {}", input_file);
            }
            Ok(proteins.len())
        })
        .collect::<Result<_>>()?;
    timings.push(("GFF processing", stage_started.elapsed()));
    log::info!(
        "Processed {} genomes, {} genes, {} threads in {:.3}s",
        all_genes.len(),
        all_genes.iter().sum::<usize>(),
        args.threads,
        stage_started.elapsed().as_secs_f64()
    );

    // Step 2: Run CD-HIT clustering
    log::info!("Step 2: Clustering sequences with CD-HIT...");
    let clustering_started = Instant::now();
    cdhit::check_cdhit_version("cd-hit")?;

    // Combine all protein sequences for clustering
    let combined_prot = outdir_path.join("combined_protein_CDS.fasta");
    let mut combined_file = std::fs::File::create(&combined_prot)?;
    let combined_dna = outdir_path.join("combined_DNA_CDS.fasta");
    let mut combined_dna_file = std::fs::File::create(&combined_dna)?;
    for i in 0..all_genes.len() {
        let prot_file = outdir_path.join(format!("input_{}.prot.fa", i));
        std::io::copy(&mut std::fs::File::open(prot_file)?, &mut combined_file)?;
        let dna_file = outdir_path.join(format!("input_{}.dna.fa", i));
        std::io::copy(&mut std::fs::File::open(dna_file)?, &mut combined_dna_file)?;
    }
    drop(combined_file);
    drop(combined_dna_file);

    let output_file = outdir_path.join("cd-hit_output.fasta");
    let data_file = outdir_path.join("gene_data.csv");

    // Run CD-HIT. Work that neither feeds cd-hit nor depends on its result
    // runs alongside it: assembling gene_data.csv (only read after
    // clustering), reading the protein-file sequence IDs, and loading the
    // annotation-ID map from the new gene_data.csv. Their results, including
    // any error, are consumed exactly where the original computed them.
    let stage_started = Instant::now();
    let (cdhit_result, early) = std::thread::scope(|scope| {
        let early = scope.spawn(|| {
            prepare_alongside_clustering(outdir_path, &data_file, all_genes.len(), &combined_prot)
        });
        let result = cdhit::run_cdhit(
            combined_prot.to_str().unwrap(),
            output_file.to_str().unwrap(),
            rustaroo_args.id.unwrap_or(0.98),
            rustaroo_args.len_dif_percent.unwrap_or(0.98),
            args.threads,
            &outdir,
        );
        (result, early.join())
    });
    let early = early.unwrap_or_else(|panic| std::panic::resume_unwind(panic));
    cdhit_result?;
    log::info!(
        "Initial clustering completed in {:.3}s",
        stage_started.elapsed().as_secs_f64()
    );
    timings.push((
        "Initial clustering and input combination",
        clustering_started.elapsed(),
    ));

    // Step 3: Generate initial network
    log::info!("Step 3: Generating pangenome network...");
    let stage_started = Instant::now();
    let cluster_file = outdir_path.join("cd-hit_output.fasta.clstr");

    // gene_data.csv was assembled alongside clustering; its errors surface here,
    // where the assembly used to run.
    let EarlyData {
        gene_data_len,
        seq_ids,
        annotation_ids,
    } = early?;
    let (annotation_lines, mut annotation_ids) = match annotation_ids {
        Ok((ids, lines)) => (lines, Some(Ok(ids))),
        Err(error) => (0, Some(Err(error))),
    };

    let prot_seq_file = combined_prot.to_str().unwrap();
    let (mut graph, centroid_context, seqid_to_centroid) = network::generate_network_with_ids(
        cluster_file.to_str().unwrap(),
        data_file.to_str().unwrap(),
        false,
        || seq_ids,
    )?;

    log::info!(
        "Initial network: {} nodes, {} edges",
        graph.node_count(),
        graph.edge_count()
    );
    timings.push(("Initial network generation", stage_started.elapsed()));

    // Step 4: Clean network
    log::info!("Step 4: Cleaning network...");

    let stage_started = Instant::now();
    clean_network::collapse_paralogs(&mut graph, &centroid_context)?;
    drop(centroid_context);
    timings.push(("Paralog context collapse", stage_started.elapsed()));
    log::info!(
        "Paralog context collapse completed in {:.3}s",
        stage_started.elapsed().as_secs_f64()
    );

    let family_started = Instant::now();
    let stage_started = Instant::now();
    let dna_index = family_index::FamilyIndex::build(
        &graph,
        &outdir,
        0.98,
        args.family_len_dif_percent,
        args.threads,
        true,
    )?;
    timings.push(("DNA similarity index", stage_started.elapsed()));
    let stage_started = Instant::now();
    clean_network::collapse_families_indexed(
        &mut graph,
        &seqid_to_centroid,
        &dna_index,
        &[1, 2, 3],
        args.quiet,
        true,
        args.length_outlier_support_proportion,
    )?;
    drop(dna_index);
    timings.push(("DNA family collapse", stage_started.elapsed()));
    let family_threshold = rustaroo_args.family_threshold.unwrap_or(0.7);
    let stage_started = Instant::now();
    let protein_index = family_index::FamilyIndex::build(
        &graph,
        &outdir,
        family_threshold,
        args.family_len_dif_percent,
        args.threads,
        false,
    )?;
    timings.push(("Protein similarity index", stage_started.elapsed()));
    let stage_started = Instant::now();
    clean_network::collapse_families_indexed(
        &mut graph,
        &seqid_to_centroid,
        &protein_index,
        &[1, 2, 3],
        args.quiet,
        false,
        args.length_outlier_support_proportion,
    )?;
    timings.push(("Protein family collapse", stage_started.elapsed()));

    // Trim trailing ends after correcting mistranslations and families.
    let stage_started = Instant::now();
    let min_support = rustaroo_args.min_trailing_support.unwrap_or(3);
    let max_recursive = rustaroo_args.trailing_recursive.unwrap_or(2);
    clean_network::trim_low_support_trailing_ends(&mut graph, min_support, max_recursive)?;
    timings.push(("Trailing-end trimming", stage_started.elapsed()));

    log::info!(
        "Family cleaning completed in {:.3}s",
        family_started.elapsed().as_secs_f64()
    );

    log::info!(
        "Cleaned network: {} nodes, {} edges",
        graph.node_count(),
        graph.edge_count()
    );

    // Step 5: Find missing genes
    let mut refind_ids = None;
    if args.refind_mode != RefindMode::Off {
        log::info!("Step 5: Finding missing genes...");
        let refind_started = Instant::now();
        let stage_started = Instant::now();
        let dna_file = outdir_path.join("combined_DNA_CDS.fasta");
        let ids = find_missing::find_missing_with_ids(
            &mut graph,
            &rustaroo_args.input_files,
            data_file.to_str().unwrap(),
            dna_file.to_str().unwrap(),
            prot_seq_file,
            &find_missing::RefindOptions {
                search_radius: args.search_radius,
                prop_match: args.refind_prop_match,
                pairwise_id_thresh: rustaroo_args.id.unwrap_or(0.98),
                merge_id_thresh: family_threshold.max(0.8),
                remove_by_consensus: rustaroo_args.remove_by_consensus.unwrap_or(false),
                only_valid_genes: args.refind_mode == RefindMode::Strict,
                creation_order_ties: args.refind_creation_order,
                verbose: !args.quiet,
            },
            annotation_ids.take(),
        )?;
        // Refinding appended to gene_data.csv; the presence/absence map is
        // this map extended with the appended lines (see below).
        refind_ids = Some(ids);
        timings.push(("Missing-gene refinding", stage_started.elapsed()));
        let stage_started = Instant::now();
        clean_network::collapse_families_indexed(
            &mut graph,
            &seqid_to_centroid,
            &protein_index,
            &[1, 2, 3],
            args.quiet,
            false,
            args.length_outlier_support_proportion,
        )?;
        timings.push(("Post-refinding family collapse", stage_started.elapsed()));
        let stage_started = Instant::now();
        clean_network::clean_misassembly_edges(
            &mut graph,
            rustaroo_args.edge_support_threshold.unwrap_or(0),
        );
        timings.push(("Misassembly edge cleaning", stage_started.elapsed()));
        log::info!(
            "Missing-gene search completed in {:.3}s",
            refind_started.elapsed().as_secs_f64()
        );
    }

    // Step 6: Generate alignments
    if !args.skip_core_aln {
        log::info!("Step 6: Generating core gene alignments...");
        let stage_started = Instant::now();
        let core_threshold = 0.99; // 99% of genomes
        alignments::generate_core_gene_alignments(
            &graph,
            &outdir,
            &args.aligner,
            args.threads,
            core_threshold,
            args.quiet,
        )?;
        timings.push(("Core gene alignment", stage_started.elapsed()));
        log::info!(
            "Core alignment completed in {:.3}s",
            stage_started.elapsed().as_secs_f64()
        );
    }

    // Step 7: Generate output files
    log::info!("Step 7: Generating output files...");
    let output_started = Instant::now();
    let stage_started = Instant::now();

    // Create isolate mapping
    let mut mems_to_isolates = std::collections::HashMap::new();
    for (i, input_file) in rustaroo_args.input_files.iter().enumerate() {
        let name = std::path::Path::new(input_file)
            .file_stem()
            .and_then(|s| s.to_str())
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("isolate_{}", i));
        mems_to_isolates.insert(i, name);
    }

    // Same as load_annotation_ids on the final gene_data.csv: the map of its
    // initial contents, extended with any lines refinding appended.
    let orig_ids = match refind_ids {
        Some(ids) => {
            find_missing::load_annotation_ids_from(
                data_file.to_str().unwrap(),
                ids,
                gene_data_len,
                annotation_lines,
            )?
            .0
        }
        None => annotation_ids
            .take()
            .expect("annotation IDs are consumed only by refinding")?,
    };

    // Generate presence/absence matrix
    let family_names = output::generate_roary_gene_presence_absence(
        &graph,
        &mems_to_isolates,
        &orig_ids,
        &outdir,
    )?;
    timings.push(("Presence/absence output", stage_started.elapsed()));

    // Generate pan-genome reference
    let stage_started = Instant::now();
    output::generate_pan_genome_reference(
        &graph,
        data_file.to_str().unwrap(),
        &rustaroo_args.input_files,
        &family_names,
        &outdir,
    )?;
    timings.push(("Reference generation", stage_started.elapsed()));

    if args.minimal {
        log::info!("Generated presence/absence matrix and reference");
    } else {
        let stage_started = Instant::now();
        let gml_file = outdir_path.join("final_graph.gml");
        output::write_gml(&graph, gml_file.to_str().unwrap())?;
        output::write_summary_stats(&graph, &outdir, all_genes.len())?;
        timings.push(("Graph and summary output", stage_started.elapsed()));
        log::info!("Generated presence/absence matrix, reference, and graph files");
    }
    log::info!(
        "Output generation completed in {:.3}s",
        output_started.elapsed().as_secs_f64()
    );

    let total = pipeline_started.elapsed();
    let measured: std::time::Duration = timings.iter().map(|(_, elapsed)| *elapsed).sum();
    timings.push(("Other pipeline overhead", total.saturating_sub(measured)));
    log::info!("Timing breakdown (wall time):");
    for (name, elapsed) in timings {
        log::info!(
            "  {:<40} {:>10.3}s ({:>5.1}%)",
            name,
            elapsed.as_secs_f64(),
            100.0 * elapsed.as_secs_f64() / total.as_secs_f64()
        );
    }
    log::info!("Pipeline total: {:.3}s", total.as_secs_f64());

    Ok(())
}

/// Results of the work run alongside the initial cd-hit clustering.
struct EarlyData {
    /// Length of gene_data.csv as first written (before refinding appends).
    gene_data_len: u64,
    seq_ids: Result<Vec<String>>,
    /// `load_annotation_ids` of the initial gene_data.csv, with its line count.
    annotation_ids: Result<(find_missing::AnnotationIds, usize)>,
}

/// Assemble gene_data.csv (as the network step used to), then read what the
/// later steps would read from it and from the protein file. Errors from
/// the assembly are returned; errors from the reads are kept for the steps
/// that would have hit them.
fn prepare_alongside_clustering(
    outdir_path: &Path,
    data_file: &Path,
    genome_count: usize,
    combined_prot: &Path,
) -> Result<EarlyData> {
    let mut data_csv = std::io::BufWriter::new(std::fs::File::create(data_file)?);
    writeln!(
        data_csv,
        "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description"
    )?;
    for i in 0..genome_count {
        let csv_file = outdir_path.join(format!("input_{}.csv", i));
        if csv_file.exists() {
            let content = std::fs::read_to_string(&csv_file)?;
            write_csv_body(&mut data_csv, &content)?;
        }
    }
    data_csv.flush()?;
    let gene_data_len = data_csv
        .into_inner()
        .map_err(|error| error.into_error())?
        .metadata()?
        .len();
    let seq_ids = network::load_sequence_ids(combined_prot.to_str().unwrap());
    let annotation_ids = find_missing::load_annotation_ids_from(
        data_file.to_str().unwrap(),
        Default::default(),
        0,
        0,
    );
    Ok(EarlyData {
        gene_data_len,
        seq_ids,
        annotation_ids,
    })
}

/// Write every line of `content` after the header, each followed by `\n`,
/// exactly as `for line in content.lines().skip(1) { writeln!(out, "{line}") }`.
/// Without any `\r`, that output is the raw text after the first newline
/// (plus a final newline if missing), so it is copied in one write.
fn write_csv_body(out: &mut impl Write, content: &str) -> io::Result<()> {
    if memchr::memchr(b'\r', content.as_bytes()).is_some() {
        for line in content.lines().skip(1) {
            writeln!(out, "{}", line)?;
        }
        return Ok(());
    }
    let Some(header_end) = memchr::memchr(b'\n', content.as_bytes()) else {
        return Ok(());
    };
    let body = &content[header_end + 1..];
    out.write_all(body.as_bytes())?;
    if !body.is_empty() && !body.ends_with('\n') {
        out.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::write_csv_body;
    use std::io::Write;

    #[test]
    fn csv_body_copy_matches_line_loop() {
        for content in [
            "", "h", "h\n", "h\n\n", "h\na", "h\na\n", "h\na\nb", "h\na\n\nb\n\n",
            "h\r\na\r\n", "h\na\rb\n", "h\na\r", "\n\n\n", "\na",
        ] {
            let mut expected = Vec::new();
            for (j, line) in content.lines().enumerate() {
                if j > 0 {
                    writeln!(expected, "{}", line).unwrap();
                }
            }
            let mut actual = Vec::new();
            write_csv_body(&mut actual, content).unwrap();
            assert_eq!(actual, expected, "{content:?}");
        }
    }
}
