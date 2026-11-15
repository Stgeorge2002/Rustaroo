//! Rustaroo: a Rust pipeline for pangenome investigation.

pub const VERSION: &str = "1.0.0";

pub mod abundance_plots;
pub mod alignment;
pub mod alignments;
pub mod args;
pub mod biocode_convert;
pub mod cdhit;
pub mod cdhit_iterative;
pub mod clean_network;
pub mod estimate_fmg;
pub mod estimate_img;
pub mod extract_gene;
pub mod family_index;
pub mod filter_pa;
pub mod find_missing;
pub mod integrate;
pub mod merge_graphs;
pub mod merge_nodes;
pub mod neighborhood;
pub mod network;
pub mod output;
pub mod post_run_gff;
pub mod post_run_msa;
pub mod prodigal;
pub mod prokka;
pub mod qc_plots;
pub mod reference;
pub mod run_prokka;
pub mod single_linkage;
pub mod spydrpick;
pub mod summary_plots;
pub mod validation;

// Re-exports for convenience
pub use args::set_default_args;
pub use validation::{is_valid_file, is_valid_folder, is_valid_gene};
