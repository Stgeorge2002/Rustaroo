//! Missing-gene sequence search in neighboring annotation windows.

use crate::alignment::{RefindMatcher, RefindScratch};
use crate::network::Graph;
use crate::network::{GeneEdge, GeneNode};
use crate::prokka::{reverse_complement, translate};
use anyhow::{Context, Result};
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::Arc;
use std::time::Instant;

/// Search for missing genes using neighboring genome membership.
pub struct RefindOptions {
    pub search_radius: usize,
    pub prop_match: f64,
    pub pairwise_id_thresh: f64,
    pub merge_id_thresh: f64,
    pub remove_by_consensus: bool,
    pub only_valid_genes: bool,
    pub verbose: bool,
}

pub fn find_missing(
    graph: &mut Graph<GeneNode, GeneEdge>,
    gff_files: &[String],
    gene_data_file: &str,
    dna_seq_file: &str,
    prot_seq_file: &str,
    options: &RefindOptions,
) -> Result<()> {
    let stage_started = Instant::now();
    if options.verbose {
        println!("Finding missing genes...");
    }

    let id_to_gff = load_annotation_ids(gene_data_file)?;
    let search_list = build_search_list(graph, &id_to_gff)?;
    let merged_sequences = load_merged_sequences(graph, gene_data_file)?;
    let preparation_elapsed = stage_started.elapsed();

    let n_searches: usize = search_list.values().map(|v| v.searches.len()).sum();
    let n_conflicts: usize = search_list.values().map(|v| v.conflicts.len()).sum();
    if options.verbose {
        println!("Number of searches to perform: {}", n_searches);
        println!("Conflict locations to resolve: {}", n_conflicts);
    }

    let progress = if options.verbose {
        let progress = ProgressBar::new((n_searches + n_conflicts) as u64);
        progress.set_style(ProgressStyle::with_template(
            "Refinding [{bar:40.cyan/blue}] {pos}/{len} tasks ({percent}%) ETA {eta_precise} {msg}",
        )?);
        progress.set_message("searching");
        progress.enable_steady_tick(std::time::Duration::from_millis(100));
        progress
    } else {
        ProgressBar::hidden()
    };

    let stage_started = Instant::now();
    let cache = RefindCache::default();
    let genome_hits: Result<Vec<_>> = gff_files
        .par_iter()
        .enumerate()
        .map(
            |(member_idx, gff_file)| match search_list.get(&member_idx) {
                Some(searches) => search_gff_file(
                    gff_file,
                    searches,
                    &merged_sequences,
                    member_idx,
                    options,
                    Some(&progress),
                    &cache,
                )
                .with_context(|| format!("Missing-gene search failed for {}", gff_file)),
                _ => Ok(GenomeHits::default()),
            },
        )
        .collect();
    let mut genome_hits = match genome_hits {
        Ok(hits) => hits,
        Err(error) => {
            progress.abandon_with_message("search failed");
            return Err(error);
        }
    };
    progress.finish_with_message("searches complete");
    drop(cache);
    let search_elapsed = stage_started.elapsed();

    let stage_started = Instant::now();
    filter_conflicts(graph, &mut genome_hits, options.remove_by_consensus)?;
    let filtering_elapsed = stage_started.elapsed();

    let stage_started = Instant::now();
    let mut n_found = 0;
    let mut dna_out = BufWriter::new(
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(dna_seq_file)?,
    );
    let mut prot_out = BufWriter::new(
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(prot_seq_file)?,
    );
    let mut data_out = BufWriter::new(
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(gene_data_file)?,
    );
    for (member_idx, result) in genome_hits.into_iter().enumerate() {
        let gff_file = &gff_files[member_idx];
        for (node_idx, hit) in result.hits {
            if hit.sequence.is_empty() {
                continue;
            }

            if graph.node_weight(node_idx).is_none() {
                continue;
            }

            let target_protein = graph[node_idx]
                .protein
                .first()
                .context("Missing refound-family protein")?;
            let hit_protein = translate_to_match(&hit.sequence, target_protein);
            let node = &mut graph[node_idx];

            if !node.members.insert(member_idx) {
                continue;
            }
            node.size = node.members.len();
            node.seq_ids
                .insert(format!("{}_refound_{}", member_idx, n_found));
            if !node.dna.contains(&hit.sequence) {
                node.dna.push(hit.sequence.clone());
            }
            if !node.protein.contains(&hit_protein) {
                node.protein.push(hit_protein.clone());
            }

            // Write to output files
            writeln!(dna_out, ">{}_refound_{}", member_idx, n_found)?;
            writeln!(dna_out, "{}", hit.sequence)?;
            writeln!(prot_out, ">{}_refound_{}", member_idx, n_found)?;
            writeln!(prot_out, "{}", hit_protein)?;

            let gff_name = std::path::Path::new(gff_file)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("");
            writeln!(
                data_out,
                "{},{},{}_refound_{},{}_refound_{},{},{},,location:{}:{};strand:{}",
                gff_name,
                hit.contig,
                member_idx,
                n_found,
                member_idx,
                n_found,
                hit_protein,
                hit.sequence,
                hit.start,
                hit.end,
                if hit.strand { '+' } else { '-' }
            )?;

            n_found += 1;
        }
    }

    dna_out.flush()?;
    prot_out.flush()?;
    data_out.flush()?;
    let update_elapsed = stage_started.elapsed();
    if options.verbose {
        println!("Number of refound genes: {}", n_found);
    }
    log::info!(
        "Refinding breakdown (wall time): preparation {:.3}s; parallel searches {:.3}s; \
         conflict/consensus filtering {:.3}s; translation, graph updates and output {:.3}s",
        preparation_elapsed.as_secs_f64(),
        search_elapsed.as_secs_f64(),
        filtering_elapsed.as_secs_f64(),
        update_elapsed.as_secs_f64()
    );

    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SearchHit {
    sequence: String,
    contig: String,
    start: usize,
    end: usize,
    strand: bool,
}

#[derive(Default)]
struct GenomeSearch {
    searches: Vec<GeneSearch>,
    conflicts: BTreeSet<(NodeIndex, String)>,
}

#[derive(Default, Debug, PartialEq, Eq)]
struct GenomeHits {
    hits: Vec<(NodeIndex, SearchHit)>,
    locations: BTreeMap<NodeIndex, SearchHit>,
}

struct GeneSearch {
    node: NodeIndex,
    sequence: String,
    anchors: BTreeSet<String>,
}

#[derive(Default)]
struct RefindCache {
    matchers: DashMap<String, Arc<RefindMatcher>>,
}

impl RefindCache {
    fn get(&self, sequence: &str) -> Option<Arc<RefindMatcher>> {
        if sequence.is_empty() {
            return None;
        }
        if let Some(matcher) = self.matchers.get(sequence) {
            return Some(Arc::clone(&matcher));
        }
        Some(Arc::clone(
            self.matchers
                .entry(sequence.to_string())
                .or_insert_with(|| {
                    Arc::new(RefindMatcher::new(sequence.as_bytes()).expect("query is nonempty"))
                })
                .value(),
        ))
    }
}

#[derive(Default)]
struct SearchScratch {
    alignment: RefindScratch,
    padded: String,
}

fn search_gff_file(
    gff_file: &str,
    search: &GenomeSearch,
    merged_sequences: &HashMap<(usize, NodeIndex), String>,
    member: usize,
    options: &RefindOptions,
    progress: Option<&ProgressBar>,
    cache: &RefindCache,
) -> Result<GenomeHits> {
    let content = std::fs::read_to_string(gff_file)?;
    let parts: Vec<&str> = content.split("##FASTA").collect();

    if parts.len() != 2 {
        anyhow::bail!("Expected one ##FASTA section when searching {}", gff_file);
    }

    // Parse FASTA sequences
    let mut contigs = HashMap::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in parts[1].lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if !current_id.is_empty() {
                contigs.insert(current_id.clone(), current_seq.clone());
            }
            current_id = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_seq = String::new();
        } else {
            current_seq.push_str(line);
        }
    }
    if !current_id.is_empty() {
        contigs.insert(current_id, current_seq);
    }

    let contigs: HashMap<_, _> = contigs
        .into_iter()
        .map(|(id, sequence)| {
            let reverse = reverse_complement(&sequence);
            (id, (sequence, reverse))
        })
        .collect();
    let mut annotations = HashMap::new();
    for line in parts[0]
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 9 {
            anyhow::bail!("Invalid GFF annotation while searching {}", gff_file);
        }
        if fields[2] != "CDS" {
            continue;
        }
        let id = fields[8]
            .split(';')
            .find_map(|attr| attr.trim().strip_prefix("ID="))
            .context("Missing CDS annotation ID in missing-gene search")?;
        let start: usize = fields[3].parse().context("Invalid GFF start coordinate")?;
        let end: usize = fields[4].parse().context("Invalid GFF end coordinate")?;
        let (sequence, _) = contigs
            .get(fields[0])
            .with_context(|| format!("Missing contig {} in {}", fields[0], gff_file))?;
        if start == 0 || end < start || end > sequence.len() {
            anyhow::bail!("Out-of-range GFF coordinates for annotation {}", id);
        }
        if annotations.insert(id, (fields[0], start, end)).is_some() {
            anyhow::bail!("Duplicate GFF annotation ID {}", id);
        }
    }
    let mut result = GenomeHits::default();

    for (_, anchor) in &search.conflicts {
        annotations
            .get(anchor.as_str())
            .with_context(|| format!("Missing conflict annotation {} in {}", anchor, gff_file))?;
    }
    let conflicts: Vec<_> = search.conflicts.iter().collect();
    let locations: Vec<_> = conflicts
        .into_par_iter()
        .map_init(SearchScratch::default, |scratch, (node, anchor)| {
            let &(contig_id, start, end) = &annotations[anchor.as_str()];
            let mut location = SearchHit {
                sequence: String::new(),
                contig: contig_id.to_string(),
                start: start - 1,
                end,
                strand: true,
            };
            if let Some(sequence) = merged_sequences.get(&(member, *node)) {
                let (contig, reverse) = &contigs[contig_id];
                let offset = start.saturating_sub(options.search_radius);
                let window_end = end.saturating_add(options.search_radius).min(contig.len());
                let min_match = (end - start) as f64 / sequence.len() as f64;
                let found = cache.get(sequence).and_then(|matcher| {
                    search_prepared(
                        sequence,
                        &matcher,
                        contig,
                        reverse,
                        offset,
                        window_end,
                        min_match,
                        options.merge_id_thresh,
                        contig_id,
                        scratch,
                    )
                });
                if let Some(found) = found {
                    location = found;
                } else {
                    location.start = offset;
                    location.end = offset;
                }
            }
            if let Some(progress) = progress {
                progress.inc(1);
            }
            (*node, location)
        })
        .collect();
    for (node, location) in locations {
        result.locations.insert(node, location);
    }

    for gene_search in &search.searches {
        for anchor in &gene_search.anchors {
            annotations.get(anchor.as_str()).with_context(|| {
                format!("Missing neighbor annotation {} in {}", anchor, gff_file)
            })?;
        }
    }

    let searched: Vec<_> = search
        .searches
        .par_iter()
        .map_init(SearchScratch::default, |scratch, gene_search| {
            let mut best_hit: Option<SearchHit> = None;
            let mut best_score = 0.0;
            let matcher = cache.get(&gene_search.sequence);
            let mut windows = HashMap::with_capacity(gene_search.anchors.len());

            for anchor in &gene_search.anchors {
                let &(contig_id, start, end) = &annotations[anchor.as_str()];
                let (contig_seq, rc_contig) = &contigs[contig_id];
                let offset = start.saturating_sub(options.search_radius);
                let window_end = end
                    .saturating_add(options.search_radius)
                    .min(contig_seq.len());
                let found = windows
                    .entry((contig_id, offset, window_end))
                    .or_insert_with(|| {
                        matcher.as_ref().and_then(|matcher| {
                            search_prepared(
                                &gene_search.sequence,
                                matcher,
                                contig_seq,
                                rc_contig,
                                offset,
                                window_end,
                                options.prop_match,
                                options.pairwise_id_thresh,
                                contig_id,
                                scratch,
                            )
                        })
                    });
                if let Some(hit) = found {
                    if hit.sequence.is_empty() {
                        continue;
                    }
                    let score = hit.sequence.len() as f64;
                    if score > best_score {
                        best_score = score;
                        best_hit = Some(hit.clone());
                    }
                }
            }

            let found = best_hit
                .filter(|hit| {
                    !options.only_valid_genes
                        || crate::validation::is_valid_gene(
                            &hit.sequence,
                            &translate_refound(&gene_search.sequence),
                        )
                })
                .map(|hit| (gene_search.node, hit));
            if let Some(progress) = progress {
                progress.inc(1);
            }
            found
        })
        .collect();

    for found in searched.into_iter().flatten() {
        result.locations.insert(found.0, found.1.clone());
        result.hits.push(found);
    }

    Ok(result)
}

#[cfg(test)]
fn search_both_strands(
    query: &str,
    contig: &str,
    reverse: &str,
    offset: usize,
    window_end: usize,
    prop_match: f64,
    threshold: f64,
    contig_id: &str,
) -> Option<SearchHit> {
    let matcher = RefindMatcher::new(query.as_bytes())?;
    search_prepared(
        query,
        &matcher,
        contig,
        reverse,
        offset,
        window_end,
        prop_match,
        threshold,
        contig_id,
        &mut SearchScratch::default(),
    )
}

fn search_prepared(
    query: &str,
    matcher: &RefindMatcher,
    contig: &str,
    reverse: &str,
    offset: usize,
    window_end: usize,
    prop_match: f64,
    threshold: f64,
    contig_id: &str,
    scratch: &mut SearchScratch,
) -> Option<SearchHit> {
    if query.is_empty() || offset >= window_end {
        return None;
    }
    static PREFIX: once_cell::sync::Lazy<regex::Regex> =
        once_cell::sync::Lazy::new(|| regex::Regex::new(r"^[ACGTEX]*N{20}").unwrap());
    static SUFFIX: once_cell::sync::Lazy<regex::Regex> =
        once_cell::sync::Lazy::new(|| regex::Regex::new(r"N{20}[ACGTEX]*$").unwrap());
    let padding = query.len() / 2;
    let targets = [
        &contig[offset..window_end],
        &reverse[contig.len() - window_end..contig.len() - offset],
    ];
    let mut best = None;
    let mut best_score = 0.0;
    for (strand, target) in targets.iter().enumerate() {
        scratch.padded.clear();
        scratch.padded.extend(std::iter::repeat('E').take(padding));
        scratch.padded.push_str(target);
        scratch.padded.extend(std::iter::repeat('E').take(padding));
        let padded = &scratch.padded;
        let Some(alignment) = matcher.align_with_scratch(padded.as_bytes(), &mut scratch.alignment)
        else {
            continue;
        };
        let mut variants = vec![Cow::Borrowed(padded.as_str())];
        if padded.contains("NNNNNNNNNNNNNNNNNNNN") {
            for pattern in [&*PREFIX, &*SUFFIX] {
                variants.push(
                    pattern.replace(padded.as_str(), |caps: &regex::Captures<'_>| {
                        "X".repeat(caps[0].len())
                    }),
                );
            }
        }
        for variant in variants {
            let variant = variant.as_ref();
            let segment = &variant[alignment.start..alignment.end];
            let n_x = segment.bytes().filter(|&b| b == b'X').count();
            let n_e = segment.bytes().filter(|&b| b == b'E').count();
            let length = segment.len() - n_x - n_e;
            let bases = segment.bytes().filter(|b| b"ACGT".contains(b)).count();
            if length as f64 / query.len() as f64 <= prop_match
                || bases as f64 / query.len() as f64 <= prop_match
            {
                continue;
            }
            let adjusted =
                alignment.distance as f64 - alignment.terminal_insertions as f64 - n_x as f64;
            let identity = 1.0 - adjusted / length as f64;
            let score = identity * length as f64;
            if identity <= threshold || score <= best_score {
                continue;
            }
            let (start, end) = if strand == 0 {
                (alignment.start, alignment.end)
            } else {
                (padded.len() - alignment.end, padded.len() - alignment.start)
            };
            let sequence = segment
                .replace(['X', 'E'], "N")
                .trim_matches('N')
                .to_string();
            best_score = score;
            best = Some(SearchHit {
                sequence,
                contig: contig_id.to_string(),
                start: offset + start.saturating_sub(padding).min(target.len()),
                end: offset + end.saturating_sub(padding).min(target.len()),
                strand: strand == 0,
            });
        }
    }
    best
}

/// Map clustering identifiers to original annotation identifiers.
pub fn load_annotation_ids(gene_data_file: &str) -> Result<HashMap<String, String>> {
    let file = File::open(gene_data_file)?;
    let reader = BufReader::new(file);
    let mut id_to_gff = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        if i == 0 {
            continue;
        }
        let line = line?;
        let fields: Vec<&str> = line.splitn(5, ',').collect();
        let id = fields.get(2).filter(|id| !id.is_empty()).with_context(|| {
            format!(
                "Missing clustering ID in {} at line {}",
                gene_data_file,
                i + 1
            )
        })?;
        let annotation = fields.get(3).filter(|id| !id.is_empty()).with_context(|| {
            format!(
                "Missing annotation ID in {} at line {}",
                gene_data_file,
                i + 1
            )
        })?;
        id_to_gff.insert(id.to_string(), annotation.to_string());
    }

    Ok(id_to_gff)
}

fn build_search_list(
    graph: &Graph<GeneNode, GeneEdge>,
    id_to_gff: &HashMap<String, String>,
) -> Result<HashMap<usize, GenomeSearch>> {
    let mut search_list: HashMap<usize, GenomeSearch> = HashMap::new();
    let mut seen = HashMap::new();

    for node_idx in graph.node_indices() {
        let node = &graph[node_idx];

        if node.dna.is_empty() {
            continue;
        }

        let query_seq = node.dna.get(node.max_len_id).with_context(|| {
            format!(
                "Invalid representative DNA index for node {}",
                node_idx.index()
            )
        })?;
        if query_seq.is_empty() {
            continue;
        }

        for neighbor_idx in graph.neighbors(node_idx) {
            let neighbor = &graph[neighbor_idx];

            for seq_id in &neighbor.seq_ids {
                let member = seq_id
                    .split('_')
                    .next()
                    .and_then(|s| s.parse::<usize>().ok())
                    .with_context(|| format!("Invalid genome ID in sequence {}", seq_id))?;

                let gff_id = id_to_gff
                    .get(seq_id)
                    .with_context(|| format!("Missing original annotation ID for {}", seq_id))?;
                let genome = search_list.entry(member).or_default();
                genome.conflicts.insert((neighbor_idx, gff_id.clone()));
                if !node.members.contains(&member) {
                    let searches = &mut genome.searches;
                    let position = *seen.entry((node_idx, member)).or_insert_with(|| {
                        let position = searches.len();
                        searches.push(GeneSearch {
                            node: node_idx,
                            sequence: query_seq.clone(),
                            anchors: BTreeSet::new(),
                        });
                        position
                    });
                    searches[position].anchors.insert(gff_id.clone());
                }
            }
        }
    }

    Ok(search_list)
}

fn load_merged_sequences(
    graph: &Graph<GeneNode, GeneEdge>,
    gene_data_file: &str,
) -> Result<HashMap<(usize, NodeIndex), String>> {
    let mut merged_ids = HashMap::new();
    for node in graph.node_indices() {
        if graph[node].centroid.len() > 1 || graph[node].merged_dna {
            for id in &graph[node].seq_ids {
                if merged_ids.insert(id.as_str(), node).is_some() {
                    anyhow::bail!("Sequence {} belongs to multiple merged nodes", id);
                }
            }
        }
    }
    if merged_ids.is_empty() {
        return Ok(HashMap::new());
    }
    let mut sequences = HashMap::new();
    let mut found = HashSet::new();
    for (line_number, line) in BufReader::new(File::open(gene_data_file)?)
        .lines()
        .enumerate()
    {
        if line_number == 0 {
            continue;
        }
        let line = line?;
        let fields: Vec<_> = line.splitn(7, ',').collect();
        if fields.len() < 6 {
            anyhow::bail!("Invalid gene metadata at line {}", line_number + 1);
        }
        let Some(&node) = merged_ids.get(fields[2]) else {
            continue;
        };
        if !found.insert(fields[2].to_string()) {
            anyhow::bail!("Duplicate gene metadata for {}", fields[2]);
        }
        let member: usize = fields[2]
            .split('_')
            .next()
            .context("Missing member in sequence ID")?
            .parse()?;
        let key = (member, node);
        if let std::collections::hash_map::Entry::Vacant(entry) = sequences.entry(key) {
            entry.insert(fields[5].to_string());
        } else {
            let representative = graph[node]
                .dna
                .get(graph[node].max_len_id)
                .context("Missing merged-node representative")?;
            sequences.insert(key, representative.clone());
        }
    }
    if found.len() != merged_ids.len() {
        anyhow::bail!("Missing gene metadata for merged nodes");
    }
    Ok(sequences)
}

#[derive(Default)]
struct Coverage {
    intervals: BTreeMap<usize, usize>,
}

impl Coverage {
    fn overlap(&self, start: usize, end: usize) -> usize {
        let mut covered = 0;
        if let Some((&left, &right)) = self.intervals.range(..start).next_back() {
            covered += right.min(end).saturating_sub(left.max(start));
        }
        for (&left, &right) in self.intervals.range(start..end) {
            covered += right.min(end) - left;
        }
        covered
    }

    fn insert(&mut self, mut start: usize, mut end: usize) {
        if start >= end {
            return;
        }
        if let Some((&left, &right)) = self.intervals.range(..=start).next_back() {
            if right >= start {
                start = left;
                end = end.max(right);
                self.intervals.remove(&left);
            }
        }
        let overlapping: Vec<_> = self
            .intervals
            .range(start..=end)
            .map(|(&left, &right)| (left, right))
            .collect();
        for (left, right) in overlapping {
            end = end.max(right);
            self.intervals.remove(&left);
        }
        self.intervals.insert(start, end);
    }
}

fn filter_conflicts(
    graph: &mut Graph<GeneNode, GeneEdge>,
    genomes: &mut [GenomeHits],
    remove_by_consensus: bool,
) -> Result<()> {
    let mut priority: Vec<_> = graph.node_indices().collect();
    priority.sort_by_key(|&node| std::cmp::Reverse((graph[node].size, node.index())));
    let rank: HashMap<_, _> = priority
        .into_iter()
        .enumerate()
        .map(|(rank, node)| (node, rank))
        .collect();

    for (member, genome) in genomes.iter_mut().enumerate() {
        let mut nodes: Vec<_> = genome.locations.keys().copied().collect();
        nodes.sort_by_key(|node| rank[node]);
        let mut coverage: HashMap<&str, Coverage> = HashMap::new();
        let mut rejected = HashSet::new();
        for node in nodes {
            let location = &genome.locations[&node];
            let max_length = *graph[node]
                .lengths
                .iter()
                .max()
                .with_context(|| format!("Missing lengths for family {}", node.index()))?;
            let occupied = coverage.entry(&location.contig).or_default();
            if occupied.overlap(location.start, location.end) * 2 >= max_length {
                if graph[node].members.contains(&member) {
                    crate::merge_nodes::remove_member_from_node(graph, node, member)?;
                }
                rejected.insert(node);
            } else {
                occupied.insert(location.start, location.end);
            }
        }
        genome.hits.retain(|(node, _)| !rejected.contains(node));
    }

    let empty: Vec<_> = graph
        .node_indices()
        .filter(|&node| graph[node].members.is_empty())
        .collect();
    for node in empty {
        crate::merge_nodes::delete_node(graph, node)?;
    }
    if remove_by_consensus {
        let mut hits = HashMap::<NodeIndex, usize>::new();
        for genome in genomes.iter() {
            for (node, _) in &genome.hits {
                if graph.node_weight(*node).is_some() {
                    *hits.entry(*node).or_default() += 1;
                }
            }
        }
        let rejected: Vec<_> = graph
            .node_indices()
            .filter(|node| hits.get(node).copied().unwrap_or(0) > graph[*node].size)
            .collect();
        for node in rejected {
            crate::merge_nodes::delete_node(graph, node)?;
        }
    }
    Ok(())
}

fn translate_to_match(dna: &str, target: &str) -> String {
    let target: HashSet<_> = target.as_bytes().windows(3).collect();
    let reverse = reverse_complement(dna);
    let mut best = String::new();
    let mut best_score = None;
    for frame in 0..3 {
        for sequence in [dna, reverse.as_str()] {
            let part = sequence.get(frame..).unwrap_or("");
            let mut padded = part.to_string();
            padded.extend(std::iter::repeat('N').take(3 - part.len() % 3));
            let protein = translate_refound(&padded);
            let triplets: HashSet<_> = protein.as_bytes().windows(3).collect();
            let score = triplets.intersection(&target).count();
            if best_score.map_or(true, |current| score > current) {
                best_score = Some(score);
                best = protein;
            }
        }
    }
    best
}

fn translate_refound(dna: &str) -> String {
    fn choices(base: &u8) -> &[u8] {
        if *base == b'N' {
            b"ACGT"
        } else {
            std::slice::from_ref(base)
        }
    }

    let dna = dna.to_ascii_uppercase();
    let mut protein = translate(&dna, 1).into_bytes();
    for (index, codon) in dna.as_bytes().chunks_exact(3).enumerate() {
        if !codon.contains(&b'N') || codon.iter().any(|base| !b"ACGTN".contains(base)) {
            continue;
        }
        let mut amino = None;
        let mut consistent = true;
        for &a in choices(&codon[0]) {
            for &b in choices(&codon[1]) {
                for &c in choices(&codon[2]) {
                    let sequence = [a, b, c];
                    let amino_acid =
                        translate(std::str::from_utf8(&sequence).expect("DNA is ASCII"), 1)
                            .as_bytes()[0];
                    if amino.is_some_and(|previous| previous != amino_acid) {
                        consistent = false;
                    }
                    amino = Some(amino_acid);
                }
            }
        }
        if consistent {
            protein[index] = amino.expect("Ambiguous codon has expansions");
        }
    }
    String::from_utf8(protein).expect("Translation is ASCII")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn query_cache_shares_prepared_matchers_across_workers() {
        let cache = RefindCache::default();
        let matchers: Vec<_> = (0..64)
            .into_par_iter()
            .map(|_| cache.get("ATGAAATAG").unwrap())
            .collect();
        assert!(matchers
            .iter()
            .all(|matcher| Arc::ptr_eq(&matchers[0], matcher)));
        assert_eq!(cache.matchers.len(), 1);
        assert!(cache.get("").is_none());
    }

    #[test]
    fn parallel_conflicts_and_searches_preserve_results_and_progress() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let gff = dir.path().join("isolate.gff");
        std::fs::write(&gff,
            "##gff-version 3\ncontig\ttest\tCDS\t1\t9\t.\t+\t0\tID=first\ncontig\ttest\tCDS\t20\t28\t.\t+\t0\tID=last\n##FASTA\n>contig\nATGAAATAGCCCCCCCCCCATGAAATAG\n")?;
        let unchanged = NodeIndex::new(8);
        let merged = NodeIndex::new(9);
        let search = GenomeSearch {
            conflicts: BTreeSet::from([
                (unchanged, "first".into()),
                (unchanged, "last".into()),
                (merged, "first".into()),
                (merged, "last".into()),
            ]),
            searches: (0..16)
                .map(|node| GeneSearch {
                    node: NodeIndex::new(20 + node),
                    sequence: if node % 2 == 0 {
                        "ATGAAATAG".into()
                    } else {
                        "CCCCCCCCC".into()
                    },
                    anchors: BTreeSet::from(["first".into(), "last".into()]),
                })
                .collect(),
        };
        let merged_sequences = HashMap::from([((0, merged), "ATGAAATAG".into())]);
        let options = RefindOptions {
            search_radius: 20,
            prop_match: 0.2,
            pairwise_id_thresh: 0.98,
            merge_id_thresh: 0.8,
            remove_by_consensus: false,
            only_valid_genes: false,
            verbose: false,
        };
        let cache = RefindCache::default();
        for strict in [false, true] {
            let strict_options = RefindOptions {
                only_valid_genes: strict,
                ..options
            };
            let mut expected = None;
            for threads in [1, 4] {
                let progress = ProgressBar::hidden();
                let result = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()?
                    .install(|| {
                        search_gff_file(
                            gff.to_str().unwrap(),
                            &search,
                            &merged_sequences,
                            0,
                            &strict_options,
                            Some(&progress),
                            &cache,
                        )
                    })?;
                assert_eq!(progress.position(), 20);
                assert_eq!(
                    (
                        result.locations[&unchanged].start,
                        result.locations[&unchanged].end
                    ),
                    (19, 28)
                );
                if let Some(expected) = &expected {
                    assert_eq!(&result, expected);
                } else {
                    expected = Some(result);
                }
            }
        }
        Ok(())
    }

    fn location(contig: &str, start: usize, end: usize) -> SearchHit {
        SearchHit {
            sequence: "A".repeat(end - start),
            contig: contig.into(),
            start,
            end,
            strand: true,
        }
    }

    #[test]
    fn coverage_counts_the_union_and_respects_boundaries() {
        let mut coverage = Coverage::default();
        coverage.insert(10, 20);
        coverage.insert(15, 25);
        coverage.insert(30, 35);
        assert_eq!(coverage.overlap(20, 32), 7);
        assert_eq!(coverage.overlap(25, 30), 0);
        coverage.insert(25, 30);
        assert_eq!(coverage.overlap(10, 35), 25);
    }

    #[test]
    fn conflict_filter_removes_lower_priority_member_and_reconnects_edges() -> Result<()> {
        let mut graph = Graph::default();
        let low = graph.add_node(GeneNode {
            size: 1,
            members: BTreeSet::from([0]),
            seq_ids: HashSet::from(["0_0_1".into()]),
            lengths: vec![10],
            ..Default::default()
        });
        let high = graph.add_node(GeneNode {
            size: 2,
            members: BTreeSet::from([1, 2]),
            lengths: vec![10],
            ..Default::default()
        });
        let left = graph.add_node(GeneNode {
            size: 1,
            members: BTreeSet::from([0]),
            ..Default::default()
        });
        let right = graph.add_node(GeneNode {
            size: 1,
            members: BTreeSet::from([0]),
            ..Default::default()
        });
        for neighbour in [left, right] {
            graph.add_edge(
                low,
                neighbour,
                GeneEdge {
                    size: 1,
                    members: BTreeSet::from([0]),
                },
            );
        }
        let mut genomes = vec![GenomeHits {
            hits: vec![(high, location("contig", 0, 10))],
            locations: BTreeMap::from([
                (low, location("contig", 0, 10)),
                (high, location("contig", 0, 10)),
            ]),
        }];
        filter_conflicts(&mut graph, &mut genomes, false)?;
        assert!(graph.node_weight(low).is_none());
        assert_eq!(genomes[0].hits.len(), 1);
        let edge = graph
            .find_edge(left, right)
            .context("Missing replacement edge")?;
        assert_eq!(graph[edge].members, BTreeSet::from([0]));
        Ok(())
    }

    #[test]
    fn consensus_removes_only_families_with_more_hits_than_members() -> Result<()> {
        let mut graph = Graph::default();
        let rejected = graph.add_node(GeneNode {
            size: 1,
            members: BTreeSet::from([2]),
            lengths: vec![10],
            ..Default::default()
        });
        let retained = graph.add_node(GeneNode {
            size: 2,
            members: BTreeSet::from([2, 3]),
            lengths: vec![10],
            ..Default::default()
        });
        let mut genomes: Vec<_> = (0..2)
            .map(|_| GenomeHits {
                hits: vec![
                    (rejected, location("first", 0, 10)),
                    (retained, location("second", 0, 10)),
                ],
                locations: BTreeMap::from([
                    (rejected, location("first", 0, 10)),
                    (retained, location("second", 0, 10)),
                ]),
            })
            .collect();
        filter_conflicts(&mut graph, &mut genomes, true)?;
        assert!(graph.node_weight(rejected).is_none());
        assert!(graph.node_weight(retained).is_some());
        Ok(())
    }

    #[test]
    fn refound_translation_uses_best_of_six_frames() {
        assert_eq!(translate_to_match("ATGAAATAG", "MK*"), "MK*X");
        assert_eq!(translate_to_match("CTATTTCAT", "MK*"), "MK*X");
        assert_eq!(translate_refound("GCNAAN"), "AX");
    }

    #[test]
    fn search_uses_panaroo_location_ties_and_wildcards() {
        for (target, query, expected, start, end, forward) in [
            ("AAAACCGGTTCCGGTT", "CCGG", "CCGG", 10, 14, true),
            ("AAAANNNNCCGGTTTT", "CCGG", "CCGG", 8, 12, true),
            ("AAAACCGGTTTT", "CCGT", "CCGGT", 4, 9, true),
        ] {
            let reverse = reverse_complement(target);
            let hit =
                search_both_strands(query, target, &reverse, 0, target.len(), 0.0, -1.0, "test")
                    .unwrap();
            assert_eq!(
                (hit.sequence.as_str(), hit.start, hit.end, hit.strand),
                (expected, start, end, forward)
            );
        }
    }

    #[test]
    fn search_scores_using_first_optimal_alignment_path() {
        let target = "ACTTGTTTACATCGGANG";
        let hit = search_both_strands(
            "CGCCGGCCN",
            target,
            &reverse_complement(target),
            0,
            target.len(),
            0.2,
            0.7,
            "test",
        )
        .unwrap();
        assert_eq!(
            (hit.sequence.as_str(), hit.start, hit.end, hit.strand),
            ("CNTCCG", 12, 18, false)
        );
    }

    #[test]
    fn search_thresholds_are_strict() {
        let target = "ACGT";
        let reverse = reverse_complement(target);
        assert!(
            search_both_strands(target, target, &reverse, 0, target.len(), 1.0, 0.9, "test")
                .is_none()
        );
        assert!(
            search_both_strands(target, target, &reverse, 0, target.len(), 0.5, 1.0, "test")
                .is_none()
        );
        assert!(
            search_both_strands(target, target, &reverse, 0, target.len(), 0.5, 0.99, "test")
                .is_some()
        );
    }

    #[test]
    fn search_list_deduplicates_same_family_and_genome() -> Result<()> {
        let mut graph = Graph::default();
        let missing = graph.add_node(GeneNode {
            dna: vec!["ATGAAATAG".into()],
            ..Default::default()
        });
        let neighbor = graph.add_node(GeneNode {
            seq_ids: HashSet::from(["0_0_1".into(), "0_0_2".into(), "1_0_1".into()]),
            ..Default::default()
        });
        graph.add_edge(
            missing,
            neighbor,
            GeneEdge {
                size: 0,
                members: Default::default(),
            },
        );
        let ids = HashMap::from([
            ("0_0_1".into(), "first".into()),
            ("0_0_2".into(), "second".into()),
            ("1_0_1".into(), "third".into()),
        ]);
        let searches = build_search_list(&graph, &ids)?;
        assert_eq!(searches[&0].searches.len(), 1);
        assert_eq!(
            searches[&0].searches[0].anchors,
            BTreeSet::from(["first".into(), "second".into()])
        );
        assert_eq!(searches[&1].searches.len(), 1);
        Ok(())
    }

    #[test]
    fn refound_gene_is_recorded_once_with_output_sequence_id() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let gff = dir.path().join("isolate.gff");
        let data = dir.path().join("gene_data.csv");
        let dna = dir.path().join("dna.fasta");
        let protein = dir.path().join("protein.fasta");
        std::fs::write(&gff, "##gff-version 3\ncontig\ttest\tCDS\t1\t3\t.\t+\t0\tID=neighbor_1\ncontig\ttest\tCDS\t13\t15\t.\t+\t0\tID=neighbor_2\n##FASTA\n>contig\nNNNATGAAATAGNNN\n")?;
        std::fs::write(&data, "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description\nisolate,contig,0_0_1,neighbor_1,MK,ATGAAATAG,,\nisolate,contig,0_0_2,neighbor_2,MK,ATGAAATAG,,\n")?;
        let mut graph = Graph::default();
        let missing = graph.add_node(GeneNode {
            dna: vec!["ATGAAATAG".into()],
            centroid: vec!["1_0_1".into()],
            protein: vec!["MK*".into()],
            members: std::collections::BTreeSet::from([1]),
            size: 1,
            lengths: vec![9],
            seq_ids: HashSet::from(["1_0_1".into()]),
            long_centroid_id: (9, "1_0_1".into()),
            ..Default::default()
        });
        let neighbor = graph.add_node(GeneNode {
            seq_ids: HashSet::from(["0_0_1".into(), "0_0_2".into()]),
            members: std::collections::BTreeSet::from([0]),
            size: 1,
            lengths: vec![3],
            ..Default::default()
        });
        graph.add_edge(
            missing,
            neighbor,
            GeneEdge {
                size: 0,
                members: Default::default(),
            },
        );
        find_missing(
            &mut graph,
            &[gff.to_str().unwrap().into()],
            data.to_str().unwrap(),
            dna.to_str().unwrap(),
            protein.to_str().unwrap(),
            &RefindOptions {
                search_radius: 5000,
                prop_match: 0.2,
                pairwise_id_thresh: 0.99,
                merge_id_thresh: 0.8,
                remove_by_consensus: false,
                only_valid_genes: false,
                verbose: false,
            },
        )?;
        assert_eq!(graph[missing].size, 2);
        assert!(graph[missing].seq_ids.contains("0_refound_0"));
        assert_eq!(std::fs::read_to_string(dna)?.matches('>').count(), 1);
        assert_eq!(
            load_annotation_ids(data.to_str().unwrap())?["0_refound_0"],
            "0_refound_0"
        );
        let index = crate::family_index::FamilyIndex {
            centroid_to_index: HashMap::from([("1_0_1".into(), 0), ("0_0_1".into(), 1)]),
            close: vec![Default::default(); 2],
        };
        let mapping = HashMap::from([
            ("1_0_1".into(), "1_0_1".into()),
            ("0_0_1".into(), "0_0_1".into()),
            ("0_0_2".into(), "0_0_1".into()),
        ]);
        crate::clean_network::collapse_families_indexed(
            &mut graph,
            &mapping,
            &index,
            &[1, 2, 3],
            true,
            false,
            0.1,
        )?;
        assert!(graph[missing].seq_ids.contains("0_refound_0"));
        Ok(())
    }

    #[test]
    fn local_search_respects_radius_preserves_all_anchors_and_reverse_coordinates() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let gff = dir.path().join("isolate.gff");
        let sequence = format!("{}CTATTTCAT{}", "C".repeat(60), "C".repeat(30));
        std::fs::write(&gff, format!("##gff-version 3\ncontig\ttest\tCDS\t1\t3\t.\t+\t0\tID=far\ncontig\ttest\tCDS\t75\t77\t.\t+\t0\tID=near\n##FASTA\n>contig\n{sequence}\n"))?;
        let search = GeneSearch {
            node: petgraph::graph::NodeIndex::new(0),
            sequence: "ATGAAATAG".into(),
            anchors: BTreeSet::from(["far".into()]),
        };
        let options = RefindOptions {
            search_radius: 20,
            prop_match: 0.8,
            pairwise_id_thresh: 0.99,
            merge_id_thresh: 0.8,
            remove_by_consensus: false,
            only_valid_genes: false,
            verbose: false,
        };
        assert!(search_gff_file(
            gff.to_str().unwrap(),
            &GenomeSearch {
                searches: vec![search],
                ..Default::default()
            },
            &HashMap::new(),
            0,
            &options,
            None,
            &RefindCache::default(),
        )?
        .hits
        .is_empty());
        let search = GeneSearch {
            node: petgraph::graph::NodeIndex::new(0),
            sequence: "ATGAAATAG".into(),
            anchors: BTreeSet::from(["far".into(), "near".into()]),
        };
        let hits = search_gff_file(
            gff.to_str().unwrap(),
            &GenomeSearch {
                searches: vec![search],
                ..Default::default()
            },
            &HashMap::new(),
            0,
            &options,
            None,
            &RefindCache::default(),
        )?
        .hits;
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].1.sequence, "ATGAAATAG");
        assert_eq!((hits[0].1.start, hits[0].1.end), (60, 69));
        assert!(!hits[0].1.strand);
        Ok(())
    }

    #[test]
    fn parallel_gene_search_keeps_search_order() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let gff = dir.path().join("isolate.gff");
        std::fs::write(
            &gff,
            "##gff-version 3\ncontig\ttest\tCDS\t1\t3\t.\t+\t0\tID=anchor\n##FASTA\n>contig\nATGAAATAG\n",
        )?;
        let anchors = BTreeSet::from(["anchor".to_string()]);
        let searches = vec![
            GeneSearch {
                node: NodeIndex::new(4),
                sequence: "ATGAAATAG".into(),
                anchors: anchors.clone(),
            },
            GeneSearch {
                node: NodeIndex::new(2),
                sequence: "ATGAAATAG".into(),
                anchors,
            },
        ];
        let progress = ProgressBar::new(2);
        progress.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        let hits = search_gff_file(
            gff.to_str().unwrap(),
            &GenomeSearch {
                searches,
                ..Default::default()
            },
            &HashMap::new(),
            0,
            &RefindOptions {
                search_radius: 5000,
                prop_match: 0.2,
                pairwise_id_thresh: 0.98,
                merge_id_thresh: 0.8,
                remove_by_consensus: false,
                only_valid_genes: false,
                verbose: false,
            },
            Some(&progress),
            &RefindCache::default(),
        )?;
        assert_eq!(
            hits.hits.iter().map(|(node, _)| *node).collect::<Vec<_>>(),
            vec![NodeIndex::new(4), NodeIndex::new(2)]
        );
        assert_eq!(progress.position(), 2);
        Ok(())
    }
}
