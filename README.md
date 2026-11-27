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

## Performance

v1.1.1 writes the same output as v1.1.0. On the six cohorts below, every file
in the output archive matched the v1.1.0 archive except `rustaroo.log`.

Times are pipeline wall time at 32 threads on Anatra `lifesci` (dual AMD EPYC
9684X), with `--clean-mode strict --remove-invalid-genes --skip-core-aln`.
Panaroo 1.8.1 and Rustaroo 1.1.0 are the published runs. Rustaroo 1.1.1 is the
same cohorts. A speedup is the baseline wall time divided by the new wall time.

| Species | Genomes | Panaroo 1.8.1 (s) | Rustaroo 1.1.0 (s) | Rustaroo 1.1.1 (s) | 1.1.0 over Panaroo | 1.1.1 over Panaroo | 1.1.1 over 1.1.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| *K. pneumoniae* | 500 | 855.38 | 319.06 | 259.9 | 2.68× | 3.29× | 1.23× |
| *K. pneumoniae* | 1,015 | 2,662.15 | 568.31 | 478.4 | 4.68× | 5.56× | 1.19× |
| *K. pneumoniae* | 3,000 | 22,884.59 | 2,284.32 | 2,095.3 | 10.02× | 10.92× | 1.09× |
| *P. aeruginosa* | 500 | 1,541.50 | 495.58 | 452.2 | 3.11× | 3.41× | 1.10× |
| *P. aeruginosa* | 1,015 | 2,293.07 | 899.45 | 704.4 | 2.55× | 3.26× | 1.28× |
| *P. aeruginosa* | 3,000 | 14,036.67 | 2,635.33 | 2,122.4 | 5.33× | 6.61× | 1.24× |

| Species | Genomes | Files compared | Identical | Differences |
| --- | ---: | ---: | ---: | ---: |
| *K. pneumoniae* | 500 | 1,512 | 1,512 | 0 |
| *K. pneumoniae* | 1,015 | 3,057 | 3,057 | 0 |
| *K. pneumoniae* | 3,000 | 9,012 | 9,012 | 0 |
| *P. aeruginosa* | 500 | 1,512 | 1,512 | 0 |
| *P. aeruginosa* | 1,015 | 3,057 | 3,057 | 0 |
| *P. aeruginosa* | 3,000 | 9,012 | 9,012 | 0 |

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
