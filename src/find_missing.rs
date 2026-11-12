//! Missing-gene sequence search in neighboring annotation windows.

use crate::alignment::{align_sequences, calculate_identity};
use crate::network::Graph;
use crate::network::{GeneEdge, GeneNode};
use crate::prokka::{reverse_complement, translate};
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// Search for missing genes using neighboring genome membership.
pub fn find_missing(
    graph: &mut Graph<GeneNode, GeneEdge>,
    gff_files: &[String],
    gene_data_file: &str,
    dna_seq_file: &str,
    prot_seq_file: &str,
    search_radius: usize,
    prop_match: f64,
    pairwise_id_thresh: f64,
    _merge_id_thresh: f64,
    verbose: bool,
) -> Result<()> {
    if verbose {
        println!("Finding missing genes...");
    }

    let id_to_gff = load_annotation_ids(gene_data_file)?;
    let search_list = build_search_list(graph, &id_to_gff)?;

    let n_searches: usize = search_list.values().map(|v| v.len()).sum();
    if verbose {
        println!("Number of searches to perform: {}", n_searches);
    }

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

    let genome_hits: Vec<_> = gff_files
        .par_iter()
        .enumerate()
        .map(
            |(member_idx, gff_file)| match search_list.get(&member_idx) {
                Some(searches) if !searches.is_empty() => search_gff_file(
                    gff_file,
                    searches,
                    search_radius,
                    prop_match,
                    pairwise_id_thresh,
                    verbose,
                )
                .with_context(|| format!("Missing-gene search failed for {}", gff_file)),
                _ => Ok(Vec::new()),
            },
        )
        .collect::<Result<_>>()?;

    for (member_idx, hits) in genome_hits.into_iter().enumerate() {
        let gff_file = &gff_files[member_idx];
        // Add found genes to graph
        for (node_idx, hit) in hits {
            if hit.sequence.is_empty() {
                continue;
            }

            if graph.node_weight(node_idx).is_none() {
                continue;
            }

            let hit_protein = translate(&hit.sequence, 11);
            let node = &mut graph[node_idx];

            if !node.members.insert(member_idx) {
                continue;
            }
            node.size = node.members.len();
            node.seq_ids
                .insert(format!("{}_refound_{}", member_idx, n_found));
            node.dna.push(hit.sequence.clone());
            node.protein.push(hit_protein.clone());

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
                "{},{},{}_refound_{},{}_refound_{},{},{},refound,location:{}:{};strand:{}",
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
    if verbose {
        println!("Number of refound genes: {}", n_found);
    }

    cleanup_empty_nodes(graph)?;
    Ok(())
}

#[derive(Debug)]
struct SearchHit {
    sequence: String,
    contig: String,
    start: usize,
    end: usize,
    strand: bool,
}

struct GeneSearch {
    node: petgraph::graph::NodeIndex,
    sequence: String,
    anchors: BTreeSet<String>,
}

fn search_gff_file(
    gff_file: &str,
    searches: &[GeneSearch],
    search_radius: usize,
    prop_match: f64,
    pairwise_id_thresh: f64,
    _verbose: bool,
) -> Result<Vec<(petgraph::graph::NodeIndex, SearchHit)>> {
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
    let mut hits = Vec::new();

    for search in searches {
        let mut best_hit: Option<SearchHit> = None;
        let mut best_score = 0.0;

        for anchor in &search.anchors {
            let &(contig_id, start, end) = annotations.get(anchor.as_str()).with_context(|| {
                format!("Missing neighbor annotation {} in {}", anchor, gff_file)
            })?;
            let (contig_seq, rc_contig) = &contigs[contig_id];
            let offset = start.saturating_sub(search_radius);
            let window_end = end.saturating_add(search_radius).min(contig_seq.len());
            let window = &contig_seq[offset..window_end];
            if let Some(hit) =
                search_in_sequence(&search.sequence, window, prop_match, pairwise_id_thresh)
            {
                let score = hit.sequence.len() as f64;
                if score > best_score {
                    best_score = score;
                    best_hit = Some(SearchHit {
                        sequence: hit.sequence,
                        contig: contig_id.to_string(),
                        start: offset + hit.start,
                        end: offset + hit.end,
                        strand: true,
                    });
                }
            }

            let reverse_window =
                &rc_contig[contig_seq.len() - window_end..contig_seq.len() - offset];
            if let Some(hit) = search_in_sequence(
                &search.sequence,
                reverse_window,
                prop_match,
                pairwise_id_thresh,
            ) {
                let score = hit.sequence.len() as f64;
                if score > best_score {
                    best_score = score;
                    best_hit = Some(SearchHit {
                        sequence: hit.sequence,
                        contig: contig_id.to_string(),
                        start: window_end - hit.end,
                        end: window_end - hit.start,
                        strand: false,
                    });
                }
            }
        }

        if let Some(hit) = best_hit {
            hits.push((search.node, hit));
        }
    }

    Ok(hits)
}

#[derive(Debug)]
struct SearchResult {
    sequence: String,
    start: usize,
    end: usize,
}

fn search_in_sequence(
    query: &str,
    target: &str,
    prop_match: f64,
    pairwise_id_thresh: f64,
) -> Option<SearchResult> {
    if query.is_empty() || target.is_empty() {
        return None;
    }

    let query_bytes = query.as_bytes();
    let target_bytes = target.as_bytes();

    // Use edlib-style HW alignment
    let max_distance = (query.len() as f64 * (1.0 - pairwise_id_thresh)).ceil() as usize;
    let result = align_sequences(query_bytes, target_bytes, max_distance);

    if result.edit_distance == usize::MAX || result.locations.is_empty() {
        return None;
    }

    let (start, end) = result.locations[0];
    let aln_length = end - start;

    // Check if alignment is long enough
    if (aln_length as f64) < (query.len() as f64 * prop_match) {
        return None;
    }

    // Check identity threshold
    let identity = calculate_identity(result.edit_distance, aln_length);
    if identity < pairwise_id_thresh {
        return None;
    }

    let found_seq = String::from_utf8_lossy(&target_bytes[start..end]).to_string();

    Some(SearchResult {
        sequence: found_seq,
        start,
        end,
    })
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
) -> Result<HashMap<usize, Vec<GeneSearch>>> {
    let mut search_list: HashMap<usize, Vec<GeneSearch>> = HashMap::new();
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

                if !node.members.contains(&member) {
                    let gff_id = id_to_gff.get(seq_id).with_context(|| {
                        format!("Missing original annotation ID for {}", seq_id)
                    })?;
                    let searches = search_list.entry(member).or_default();
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

fn cleanup_empty_nodes(graph: &mut Graph<GeneNode, GeneEdge>) -> Result<()> {
    let mut bad_nodes = Vec::new();

    for node_idx in graph.node_indices() {
        if graph[node_idx].members.is_empty() {
            bad_nodes.push(node_idx);
        }
    }

    for node in bad_nodes {
        graph.remove_node(node);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_search_in_sequence() {
        let query = "ATGGCCTAG";
        let target = "NNNATGGCCTAGNNNN";
        let result = search_in_sequence(query, target, 0.8, 0.95);
        assert!(result.is_some());
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
        assert_eq!(searches[&0].len(), 1);
        assert_eq!(
            searches[&0][0].anchors,
            BTreeSet::from(["first".into(), "second".into()])
        );
        assert_eq!(searches[&1].len(), 1);
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
            seq_ids: HashSet::from(["1_0_1".into()]),
            ..Default::default()
        });
        let neighbor = graph.add_node(GeneNode {
            seq_ids: HashSet::from(["0_0_1".into(), "0_0_2".into()]),
            members: std::collections::BTreeSet::from([0]),
            size: 1,
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
            5000,
            0.2,
            0.99,
            0.7,
            false,
        )?;
        assert_eq!(graph[missing].size, 2);
        assert!(graph[missing].seq_ids.contains("0_refound_0"));
        assert_eq!(std::fs::read_to_string(dna)?.matches('>').count(), 1);
        assert_eq!(
            load_annotation_ids(data.to_str().unwrap())?["0_refound_0"],
            "0_refound_0"
        );
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
        assert!(
            search_gff_file(gff.to_str().unwrap(), &[search], 20, 0.8, 0.99, false)?.is_empty()
        );
        let search = GeneSearch {
            node: petgraph::graph::NodeIndex::new(0),
            sequence: "ATGAAATAG".into(),
            anchors: BTreeSet::from(["far".into(), "near".into()]),
        };
        let hits = search_gff_file(gff.to_str().unwrap(), &[search], 20, 0.8, 0.99, false)?;
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].1.sequence, "ATGAAATAG");
        assert_eq!((hits[0].1.start, hits[0].1.end), (60, 69));
        assert!(!hits[0].1.strand);
        Ok(())
    }
}
