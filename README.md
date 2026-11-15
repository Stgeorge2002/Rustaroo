# Rustaroo

Rustaroo is a Rust pipeline for pangenome investigation. The core pipeline is runnable, but this is still a work in progress: exact behavioral parity with the baseline implementation, all original utilities, and full-pipeline performance equivalence have not been established. It incorporates reference-selection code from the earlier [rustaroo-alpha](https://github.com/Stgeorge2002/rustaroo-alpha) project.

## Overview

The implementation uses Rust for annotation processing and graph operations and invokes external tools such as CD-HIT and MAFFT. The project builds the `rustaroo` executable and exposes a `rustaroo` Rust library.

## Project Structure

```
rustaroo/
├── Cargo.toml                  # Package manifest
├── Cargo.lock                  # Reproducible dependency versions
├── src/                        # Library, CLI, and reference selection
└── tests/                      # Integration and regression tests
```

## Features

### Core Functionality

- **Network Generation**: Build pangenome networks from gene annotations
- **Graph Cleaning**: Remove contamination and collapse gene families
- **Missing Gene Search**: Find genes missed in initial annotation
- **MSA Generation**: Create multiple sequence alignments
- **Summary Statistics**: Report family and graph counts

Ancillary modules, missing-gene search, and alignment generation are not yet validated as faithful replacements for the corresponding baseline workflows.

### Modes

- **Strict**: High confidence, removes likely contaminants
- **Moderate**: Balanced approach for most datasets
- **Sensitive**: Keep rare genes, useful for plasmid analysis

## Building

### Prerequisites

- Rust stable (install from [rustup.rs](https://rustup.rs))
- External tools:
  - CD-HIT (for sequence clustering)
  - MAFFT or PRANK (for alignments)
  - Optional: mash, prodigal, prokka

### Build Commands

```bash
# Build the project
cargo build --release

# Run core unit and graph regression tests
cargo test --lib --test test_network --test test_cleaning_reference
```

## Usage

### Container package

The publishing workflow builds a Linux x86-64 container containing Rustaroo,
CD-HIT and MAFFT. After the package has been published, run it with Docker:

```bash
docker run --rm ghcr.io/stgeorge2002/rustaroo:1.0.0 --version
docker run --rm -v "$PWD:/data" ghcr.io/stgeorge2002/rustaroo:1.0.0 \
  -i /data/genome1.gff /data/genome2.gff \
  -o /data/results --clean-mode strict --threads 20 --skip-core-aln
```

For an HPC using Apptainer:

```bash
apptainer pull rustaroo-1.0.0.sif docker://ghcr.io/stgeorge2002/rustaroo:1.0.0
apptainer exec --cleanenv --bind "$PWD:/data" rustaroo-1.0.0.sif \
  rustaroo -i /data/genome1.gff /data/genome2.gff \
  -o /data/results --clean-mode strict --threads 20 --skip-core-aln
```

Run analyses inside a scheduler allocation. The image does not include optional
annotation tools such as Prokka or Prodigal.

Maintainers can publish an existing release using **Actions > Publish container >
Run workflow**, selecting the default branch and entering its tag (for example
`v1.0.0`). The workflow uses packaging files from the selected workflow revision
and pipeline source from that release tag, verifies the version, and smoke-tests
the image before publishing. Future stable `vX.Y.Z` tag pushes also trigger it.
No personal access token is required: the workflow uses `GITHUB_TOKEN` with
package-write permission. If the new GHCR package is private, change its package
visibility to public to allow unauthenticated pulls. Record the published image
digest when pinning a benchmark environment.

### Main Pipeline

```bash
# Basic usage
./target/release/rustaroo -i *.gff -o output_dir --clean-mode strict

# With custom parameters
./target/release/rustaroo \
  -i genome1.gff genome2.gff genome3.gff \
  -o results \
  --clean-mode moderate \
  --threads 8 \
  -c 0.98 \
  -f 0.7 \
  --aligner mafft

# Sensitive mode for plasmid analysis
./target/release/rustaroo -i *.gff -o output --clean-mode sensitive
```

### Controlled 50-genome cleaning run

From this directory, with the input folder in its parent:

```bash
cargo build --release
./target/release/rustaroo \
  -i ../test_gff_50/*.gff3 \
  -o ../rustaroo_50_output \
  --clean-mode moderate \
  --threads 20 \
  --refind-mode off \
  --skip-core-aln
```

Use a fresh output directory. These two final flags explicitly exclude refinding and alignment; they are not equivalent to enabling those stages. Rustaroo currently requests core alignments by default. Use `./target/release/rustaroo --help` for the actual supported CLI.

To enable missing-gene search while retaining the validated no-alignment setting, omit `--refind-mode off`:

```bash
./target/release/rustaroo \
  -i ../test_gff_50/*.gff3 \
  -o ../rustaroo_50_refind_output \
  --clean-mode moderate --threads 20 --skip-core-aln
```

## External Dependencies

External tools:

- **Required**:
  - CD-HIT (v4.6+): Sequence clustering
  
- **Recommended**:
  - MAFFT or PRANK: Multiple sequence alignment
  - Prokka: Gene annotation (if starting from assemblies)

- **Optional**:
  - mash: Distance calculations
  - prodigal: Alternative gene prediction

Install on Ubuntu/Debian:
```bash
sudo apt install cd-hit mafft prokka
```

Install via conda:
```bash
conda install -c bioconda cd-hit mafft prokka
```

## Output Files

The main pipeline writes these files; exact format/metadata parity is not guaranteed:

- `gene_presence_absence.csv`: PA matrix
- `gene_data.csv`: Gene metadata
- `gene_presence_absence_roary.csv`: Extended PA matrix
- `gene_presence_absence.Rtab`: Binary PA matrix
- `pan_genome_reference.fa`: Representative sequences
- `reference_report.tsv`: Reference selection reasons, sample support, and source sequence IDs
- `final_graph.gml`: Pangenome graph
- `summary_statistics.txt`: Family and graph counts

Isolate columns follow input order. Presence/absence sequence identifiers are mapped back to annotation identifiers from `gene_data.csv`.

The GML writer includes node labels, centroid sequences, sequence/genome memberships, lengths, and edge genome memberships. Lists use NetworkX's GML encoding (including singleton markers), and strings use GML character entities rather than backslash quote escapes. Remaining graph/output parity limitations are listed below.

### Rustaroo reference selection

The final reference FASTA now uses the reference-selection policy from
[rustaroo-alpha](https://github.com/Stgeorge2002/rustaroo-alpha), pinned to revision
`6b09ae5a5a401003b6a24195132c15524109e965`. The original sequence selector,
genetic-code helpers, sequence fitting, and reference evidence rules are
ported with their MIT notice in the reference module.
Only optional profiling/hash-seed instrumentation and module paths were adapted.

This is an output-only stage after the existing pipeline. It receives
an immutable final graph: clustering, cleaning, refinding, paralog membership,
edges, and alignment inputs are not changed. Rustaroo's graph construction,
family reconciliation, quality filtering, and refinding are not used.

Selection considers every accepted member sequence in each final family,
not just its centroids:

- Count distinct samples for each normalized observed allele.
- Prefer the longest acceptable intact coding sequence, with Rustaroo's
  terminal-extension, internal-insertion, and fused-tail checks.
- Use supported, unclipped annotated gene ends to reject private extensions.
  Refound copies and multi-fragment loci do not establish finished ends.
- Exclude refound-only candidates below 95% of the longest intact annotated
  allele, without removing their memberships.
- Break equal-length/frequency ties using an observed haplotype's distance
  to the weighted codon majority; the majority itself is not emitted.
- Apply Rustaroo's read-only physical-locus tests and guarded fragment joins
  within existing families only. Synthetic models never become graph members.
- Retain a flagged `no_intact_copy` reference when no intact candidate exists.

GFF coordinates, strand, and contig lengths are read separately as reference
evidence; they do not replace the existing input parser or reinstate rejected
annotations. Refound zero-based half-open coordinates are converted to
one-based inclusive coordinates for Rustaroo's evidence rules. Reference
translation uses Rustaroo's bacterial code-11 rules; upstream translation is
unchanged.

FASTA headers reuse the exact family IDs from the presence/absence writer.
There is one record per final family, including unnamed families and separate
paralogs that share a centroid. DNA is wrapped at 60 bases, as in Rustaroo.
The former first-centroid reference writer and centroid-based skipping have
been removed. These FASTA sequences, identifiers, and record counts intentionally
differ from older versions' reference outputs.

`reference_report.tsv` records integrity, selection reason, finished-end
support, exact-reference sample counts, distinct/intact allele counts, and
source member IDs (the parts for selected synthetic joins). Missing or
ambiguous member evidence causes an explicit error, not a centroid fallback.
No family is removed because it has a weak reference.

“Same reference policy” does not mean identical outputs from running the two
whole programs: their final families and accepted member evidence differ.
Rustaroo's earlier family-repair provenance is not imported; physical-locus
evidence is recomputed read-only over the finalized members.

Validation on the 50-GFF fixture produced 2,916 references. Every sequence
and selection decision matched a harness running Rustaroo's unchanged original
selector on the same finalized-family evidence. Biological family
memberships and edge payloads matched the saved pre-change run after normalizing
existing run-dependent refound identifiers. A separate finalized-graph test
verified identical node/edge payloads, GML bytes, and all three presence/absence
files before and after reference generation. The release-mode reference stage
took approximately 0.94 seconds in that test; this is a single local measurement,
not a general performance guarantee.

## Performance

### Cleaning architecture

- Process independent GFF files in parallel, preserving input and genomic order.
- Use a stable, undirected graph so deletion does not relocate other live nodes.
- Preserve centroid/DNA/protein associations during merges, including identical sequences belonging to distinct centroids.
- Collapse paralogs using genomic context.
- Globally precluster centroid DNA, then protein, using the reference CD-HIT cascades.
- Score only candidate pairs inside preclusters, in parallel with the requested thread count.
- Reuse each query's bit-parallel Myers matcher and reverse complement.
- Store accepted similarities sparsely and reuse them at neighborhood depths 1, 2, and 3.
- Project sparse similarities into neighborhoods using union-find, retaining genome-overlap and split-gene conflict checks.

There are no CD-HIT subprocesses or sequence alignments inside the neighborhood-cleaning loop.

### Measured scope and output differences

A single 20-thread comparison with the baseline implementation on the 50-GFF dataset used the **identical** starting graph (5,104 nodes, 6,625 edges) for both implementations:

| Measurement | Result |
|---|---:|
| Standard DNA + protein cleaning | 36.320 s |
| Rust DNA + protein cleaning | 8.448 s |
| Cleaning-only speedup | 4.30x |
| Reference families not exactly reproduced | 35 / 3,092 (1.13%) |
| Genes changing family after optimal one-to-one family relabeling | 88 / 107,067 (0.082%) |
| Reference genome-labelled edge records not exactly reproduced | 142 / 3,599 (3.95%) |

Gene-assignment disagreement counts a gene as unchanged only when its reference family is matched one-to-one to a native family containing that gene; missing genes also count as disagreements. This is not a raw node-count comparison. Edge records include complete endpoint family memberships, so changed families can affect several edge records.

These single-run measurements are not a guarantee for other datasets. All 107,067 extracted input records matched standard annotation IDs, proteins, and DNA. Different disagreement metrics describe different properties; a small gene-assignment percentage does not imply an equally small difference in complete families or edge records.

Final 20-thread CLI runs also compared annotation processing, clustering, cleaning, refinding, and output generation, with alignments disabled in both implementations:

| Measurement | Result |
|---|---:|
| Standard pipeline, refinding enabled | 59.381 s |
| Rust pipeline, refinding enabled (final isolated repeat) | 14.662 s |
| Refinding-enabled pipeline speedup | 4.05x |
| Final gene-assignment disagreements, refinding enabled | 1,565 / 108,755 (1.439%) |
| Reference families not exactly reproduced, refinding enabled | 199 / 2,911 (6.84%) |
| Final gene-assignment disagreements, refinding disabled | 110 / 106,893 (0.103%) |
| Reference families not exactly reproduced, refinding disabled | 40 / 2,918 (1.37%) |

CLI durations above use a monotonic clock. With refinding disabled, both final outputs retained exactly the same annotated genes. For the refinding-enabled comparison, generated identifiers were normalized by isolate, contig, start/end coordinates, and DNA sequence, including duplicate occurrences. Memberships were read from the Roary output to avoid treating reference QC suffixes as different genes.

These are diagnostic measurements, not an approved fixed error budget. Refound calls, QC annotations, graph metadata, and alignment outputs are not claimed equivalent. No comparison with core alignment enabled has been made. Final exact-match/GML improvements and a subsequent isolated repeat preserved normalized family memberships. NetworkX successfully read the final GML, whose 107,796 gene IDs exactly matched the presence/absence output.

The timing table predates Rustaroo-based reference selection; it does not
include the additional final allele/evidence selection work.

### Missing-gene alignment

The sequence search now scans targets with bit-parallel Myers and recovers the start in a window of at most `query_length + edit_distance`, rather than allocating whole-contig quadratic DP/traceback tables. For query length `m` and target length `n`, this takes approximately `O(n * ceil(m / 64) + m * (m + d))` time and `O(m + d)` auxiliary memory, where `d` is the accepted edit distance. Regression tests compare distances, start/end coordinates, and tie-breaking against the former full-table implementation.

For an exact match, the start is simply `end - query_length`; the recovery DP and its allocations are skipped entirely without changing coordinates or tie-breaking.

Duplicate family/genome searches are grouped while preserving every distinct neighbor annotation anchor. Genomes are searched in parallel using the CLI's bounded Rayon pool; results are applied in input order. Reverse complements are cached per contig.

Searches use the reference's neighbor-local window, clipped to the contig boundaries using `--search-radius`, rather than scanning unrelated contigs. This is an intentional correction to the former whole-contig behavior. In the 50-GFF runs, the search stage dropped from 28.219 s for whole-contig scanning to 0.491 s for local windows; the refound calls changed as expected. Regression tests cover radius boundaries, multiple anchors, and reverse-strand sequence/coordinate consistency.

Refound sequence IDs are recorded in the graph and each newly found genome membership is counted once.

## Development Status

### Known limitations

- Family merging can differ due to traversal and tie ordering.
- Refinding uses local windows but reference overlap/conflict filtering, partial/ambiguous hits, translation-frame selection, consensus removal, and strict-mode semantics are incomplete.
- Non-default genetic codes, edge filtering, optional final paralog merging, GML metadata, and alignment outputs have not achieved full reference parity.
- The complete integration suite has pre-existing missing XML helper exports in validation tests. The targeted core suite below does not depend on those exports.
- No universal speedup, memory reduction, or exact output-equivalence claim is made.

## Testing

```bash
# Run all tests
cargo test

# Run validated core selectors together
cargo test --lib --test test_network --test test_cleaning_reference

# Run with output
cargo test -- --nocapture

```

The ignored reference tests use `RUSTAROO_REFERENCE_FIXTURE`, or
`RUSTAROO_REFERENCE_INPUT` and `RUSTAROO_REFERENCE_OUTPUT` with optional
`RUSTAROO_REFERENCE_THREADS`.

## Contributing

Contributions and reproducible validation results are welcome.

## Citation

For the foundational pangenome pipeline research, see:

> Tonkin-Hill, G., MacAlasdair, N., Ruis, C. et al. Genome Biol 21, 180 (2020). https://doi.org/10.1186/s13059-020-02090-4

## License

MIT License; see [LICENSE](LICENSE) for the original project notice. The
reference module also includes the rustaroo-alpha reference-selection code notice.

## Authors

- See the license notices and citation above.

## Links

- Earlier reference-selection project: https://github.com/Stgeorge2002/rustaroo-alpha
- Foundational research: https://doi.org/10.1186/s13059-020-02090-4
