# Rustaroo

Rustaroo is a Rust pipeline for pangenome analysis. It processes GFF annotations,
builds and cleans a gene-family graph, searches for missing genes, selects family
reference sequences, and can generate alignments.

## Build and run

Install stable Rust and CD-HIT. Install MAFFT or PRANK if you want alignments.

```bash
cargo build --release
cargo test --lib --test test_network --test test_cleaning_reference

./target/release/rustaroo \
  -i genome1.gff genome2.gff \
  -o results \
  --clean-mode strict \
  --threads 8 \
  --skip-core-aln
```

Cleaning modes are `strict`, `moderate`, and `sensitive`. Core alignments are
enabled by default; use `--skip-core-aln` to omit them.

## Results

The output directory includes presence/absence tables, family references,
`reference_report.tsv`, `final_graph.gml`, and summary statistics. Use
`--minimal` to retain only the log, reference FASTA, and main presence/absence
table. Output and reference behavior may differ from other implementations.

## Tests

Run the full suite with `cargo test`. Reference tests may require fixtures set
through `RUSTAROO_REFERENCE_FIXTURE` or `RUSTAROO_REFERENCE_INPUT` and
`RUSTAROO_REFERENCE_OUTPUT`.

## Citation and license

Tonkin-Hill et al. *Genome Biology* 21, 180 (2020).
<https://doi.org/10.1186/s13059-020-02090-4>

MIT License; see [LICENSE](LICENSE). Reference-selection code includes its own
license notice.
