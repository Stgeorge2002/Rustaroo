//! Generate output files - presence/absence matrices, GML, alignments
//!
//! Port of core output generation from rustaroo/generate_output.py

use crate::network::Graph;
use crate::network::{GeneEdge, GeneNode};
use anyhow::{Context, Result};
use petgraph::graph::NodeIndex;
use petgraph::visit::{IntoEdgeReferences, IntoNodeIdentifiers};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub use crate::reference::{
    generate_pan_genome_reference, generate_pan_genome_reference_with, preload_reference,
};

/// Generate Roary-style gene presence/absence CSV
pub fn generate_roary_gene_presence_absence(
    graph: &Graph<GeneNode, GeneEdge>,
    mems_to_isolates: &HashMap<usize, String>,
    orig_ids: &HashMap<String, String, impl std::hash::BuildHasher>,
    output_dir: &str,
) -> Result<HashMap<NodeIndex, String>> {
    let mut members: Vec<_> = mems_to_isolates.iter().collect();
    members.sort_by_key(|(member, _)| **member);
    let isolates: Vec<_> = members.iter().map(|(_, name)| (*name).clone()).collect();
    let mems_to_index: HashMap<usize, usize> = members
        .into_iter()
        .enumerate()
        .map(|(i, (mem, _))| (*mem, i))
        .collect();

    // Create output files
    let roary_path = Path::new(output_dir).join("gene_presence_absence_roary.csv");
    let csv_path = Path::new(output_dir).join("gene_presence_absence.csv");
    let rtab_path = Path::new(output_dir).join("gene_presence_absence.Rtab");

    let mut roary_file = BufWriter::new(File::create(roary_path)?);
    let mut csv_file = BufWriter::new(File::create(csv_path)?);
    let mut rtab_file = BufWriter::new(File::create(rtab_path)?);

    // Write headers
    let mut header = vec![
        "Gene".to_string(),
        "Non-unique Gene name".to_string(),
        "Annotation".to_string(),
        "No. isolates".to_string(),
        "No. sequences".to_string(),
        "Avg sequences per isolate".to_string(),
        "Genome Fragment".to_string(),
        "Order within Fragment".to_string(),
        "Accessory Fragment".to_string(),
        "Accessory Order with Fragment".to_string(),
        "QC".to_string(),
        "Min group size nuc".to_string(),
        "Max group size nuc".to_string(),
        "Avg group size nuc".to_string(),
    ];
    header.extend(isolates.clone());

    writeln!(roary_file, "{}", header.join(","))?;

    let simple_header = vec!["Gene", "Non-unique Gene name", "Annotation"];
    let mut csv_header = simple_header
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();
    csv_header.extend(isolates.clone());
    writeln!(csv_file, "{}", csv_header.join(","))?;

    let mut rtab_header = vec!["Gene".to_string()];
    rtab_header.extend(isolates.clone());
    writeln!(rtab_file, "{}", rtab_header.join("\t"))?;

    // Process each connected component (gene family)
    let mut used_gene_names = HashSet::new();
    used_gene_names.insert("".to_string());
    let mut unique_id_count = 0;
    let mut frag = 0;
    let mut entries = Vec::new();
    let mut family_names = HashMap::new();

    // Simple connected component detection
    for node_idx in graph.node_identifiers() {
        frag += 1;
        let node = &graph[node_idx];

        let name = node
            .annotation
            .trim()
            .trim_matches(';')
            .split(';')
            .filter(|s| !s.is_empty())
            .collect::<Vec<_>>()
            .join("~~~");

        let clean_name: String = name
            .chars()
            .filter(|c| c.is_alphanumeric() || *c == '_' || *c == '~')
            .collect();

        let gene_name = if !used_gene_names.contains(&clean_name.to_lowercase()) {
            used_gene_names.insert(clean_name.to_lowercase());
            clean_name
        } else {
            let fallback = format!("group_{}", unique_id_count);
            unique_id_count += 1;
            fallback
        };
        family_names.insert(node_idx, gene_name.clone());

        // Build presence/absence vector
        let mut pres_abs = BTreeMap::<usize, String>::new();
        let mut seq_ids: Vec<_> = node.seq_ids.iter().collect();
        seq_ids.sort();
        for seq_id in seq_ids {
            let (member, _) = seq_id.split_once('_').with_context(|| {
                format!(
                    "Invalid clustering ID in presence/absence output: {}",
                    seq_id
                )
            })?;
            let member: usize = member
                .parse()
                .with_context(|| format!("Invalid genome ID in sequence {}", seq_id))?;
            let idx = *mems_to_index
                .get(&member)
                .with_context(|| format!("Missing isolate name for genome {}", member))?;
            let id = orig_ids.get(seq_id).unwrap_or(seq_id);
            let value = pres_abs.entry(idx).or_default();
            if value.is_empty() {
                *value = id.clone();
            } else {
                value.push_str(&format!(";{}", id));
            }
        }

        let entry_size = pres_abs.values().filter(|s| !s.is_empty()).count();

        let entry = vec![
            gene_name.clone(),
            node.annotation.clone(),
            node.description.clone(),
            node.size.to_string(),
            node.seq_ids.len().to_string(),
            format!("{:.2}", node.seq_ids.len() as f64 / node.size.max(1) as f64),
            frag.to_string(),
            "1".to_string(),
            "".to_string(),
            "".to_string(),
            "".to_string(),
            node.lengths.iter().min().unwrap_or(&0).to_string(),
            node.lengths.iter().max().unwrap_or(&0).to_string(),
            format!(
                "{:.1}",
                if !node.lengths.is_empty() {
                    node.lengths.iter().sum::<usize>() as f64 / node.lengths.len() as f64
                } else {
                    0.0
                }
            ),
        ];
        entries.push((entry_size, gene_name, entry, pres_abs));
    }

    // Sort by entry size (most common first)
    entries.sort_by(|a, b| b.0.cmp(&a.0));

    for (_size, gene_name, mut entry, sparse) in entries {
        let mut pres_abs = vec!["".to_string(); isolates.len()];
        for (index, value) in sparse {
            pres_abs[index] = value;
        }
        entry.extend(pres_abs.iter().cloned());
        writeln!(roary_file, "{}", entry.join(","))?;

        let mut simple_entry = entry[0..3].to_vec();
        simple_entry.extend(pres_abs.iter().map(|s| s.to_string()));
        writeln!(csv_file, "{}", simple_entry.join(","))?;

        write!(rtab_file, "{}\t", gene_name)?;
        let binary: Vec<_> = pres_abs
            .iter()
            .map(|s| if s.is_empty() { "0" } else { "1" })
            .collect();
        writeln!(rtab_file, "{}", binary.join("\t"))?;
    }

    roary_file.flush()?;
    csv_file.flush()?;
    rtab_file.flush()?;
    Ok(family_names)
}

/// Write GML graph file (NetworkX-compatible format)
pub fn write_gml(graph: &Graph<GeneNode, GeneEdge>, output_path: &str) -> Result<()> {
    let mut file = BufWriter::new(File::create(output_path)?);

    writeln!(file, "graph [")?;
    writeln!(file, "  directed 0")?;

    // Write nodes
    for node_idx in graph.node_identifiers() {
        let node = &graph[node_idx];
        writeln!(file, "  node [")?;
        writeln!(file, "    id {}", node_idx.index())?;
        writeln!(file, "    label \"{}\"", node_idx.index())?;
        writeln!(file, "    size {}", node.size)?;
        writeln!(file, "    maxLenId {}", node.max_len_id)?;
        writeln!(file, "    hasEnd {}", if node.has_end { 1 } else { 0 })?;
        writeln!(file, "    paralog {}", if node.paralog { 1 } else { 0 })?;
        writeln!(
            file,
            "    mergedDNA {}",
            if node.merged_dna { 1 } else { 0 }
        )?;
        write_gml_list(&mut file, "centroid", &node.centroid, true)?;
        write_gml_list(
            &mut file,
            "members",
            &node.members.iter().collect::<Vec<_>>(),
            false,
        )?;
        let mut seq_ids: Vec<_> = node.seq_ids.iter().collect();
        seq_ids.sort();
        write_gml_list(&mut file, "seqIDs", &seq_ids, true)?;
        write_gml_list(&mut file, "dna", &node.dna, true)?;
        write_gml_list(&mut file, "protein", &node.protein, true)?;
        write_gml_list(&mut file, "lengths", &node.lengths, false)?;
        writeln!(file, "    longCentroidID {}", node.long_centroid_id.0)?;
        writeln!(
            file,
            "    longCentroidID \"{}\"",
            escape_gml(&node.long_centroid_id.1)
        )?;
        writeln!(file, "    annotation \"{}\"", escape_gml(&node.annotation))?;
        writeln!(
            file,
            "    description \"{}\"",
            escape_gml(&node.description)
        )?;
        writeln!(file, "  ]")?;
    }

    // Write edges
    use petgraph::visit::EdgeRef;
    for edge in graph.edge_references() {
        writeln!(file, "  edge [")?;
        writeln!(file, "    source {}", edge.source().index())?;
        writeln!(file, "    target {}", edge.target().index())?;
        writeln!(file, "    size {}", edge.weight().size)?;
        write_gml_list(
            &mut file,
            "members",
            &edge.weight().members.iter().collect::<Vec<_>>(),
            false,
        )?;
        writeln!(file, "  ]")?;
    }

    writeln!(file, "]")?;
    file.flush()?;
    Ok(())
}

/// Escape special characters for GML format
fn escape_gml(s: &str) -> String {
    let mut escaped = String::new();
    for character in s.chars() {
        if matches!(character, '&' | '"' | '\r' | '\n') || !character.is_ascii() {
            escaped.push_str(&format!("&#{};", character as u32));
        } else {
            escaped.push(character);
        }
    }
    escaped
}

fn write_gml_list<T: std::fmt::Display>(
    file: &mut impl Write,
    key: &str,
    values: &[T],
    quoted: bool,
) -> Result<()> {
    match values.len() {
        0 => writeln!(file, "    {} \"[]\"", key)?,
        1 => writeln!(file, "    {} \"_networkx_list_start\"", key)?,
        _ => {}
    }
    for value in values {
        if quoted {
            writeln!(file, "    {} \"{}\"", key, escape_gml(&value.to_string()))?;
        } else {
            writeln!(file, "    {} {}", key, value)?;
        }
    }
    Ok(())
}

/// Generate summary statistics
pub fn write_summary_stats(
    graph: &Graph<GeneNode, GeneEdge>,
    output_dir: &str,
    n_isolates: usize,
) -> Result<()> {
    let stats_path = Path::new(output_dir).join("summary_statistics.txt");
    let mut file = BufWriter::new(File::create(stats_path)?);

    let n_nodes = graph.node_count();
    let n_edges = graph.edge_count();

    let core_genes = graph
        .node_identifiers()
        .filter(|&n| graph[n].size >= n_isolates)
        .count();

    let accessory_genes = n_nodes - core_genes;

    writeln!(file, "Number of isolates: {}", n_isolates)?;
    writeln!(file, "Total gene families: {}", n_nodes)?;
    writeln!(file, "Core gene families (>= 99%): {}", core_genes)?;
    writeln!(file, "Accessory gene families: {}", accessory_genes)?;
    writeln!(file, "Total edges in graph: {}", n_edges)?;
    file.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_escape_gml() {
        assert_eq!(escape_gml("test\"quote"), "test&#34;quote");
        assert_eq!(escape_gml("test\\slash"), "test\\slash");
        assert_eq!(escape_gml("&\n\u{03b1}"), "&#38;&#10;&#945;");
    }

    #[test]
    fn gml_preserves_family_membership_sequences_and_singleton_lists() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let path = dir.path().join("graph.gml");
        let mut graph = Graph::default();
        let first = graph.add_node(GeneNode {
            size: 1,
            centroid: vec!["0_0_1".into()],
            members: std::collections::BTreeSet::from([0]),
            seq_ids: HashSet::from(["0_0_1".into()]),
            dna: vec!["ATGAAA".into()],
            protein: vec!["MK".into()],
            lengths: vec![6],
            long_centroid_id: (6, "0_0_1".into()),
            annotation: "a\"b&c".into(),
            ..Default::default()
        });
        let second = graph.add_node(graph[first].clone());
        graph.add_edge(
            first,
            second,
            GeneEdge {
                size: 1,
                members: std::collections::BTreeSet::from([0]),
                ..GeneEdge::default()
            },
        );
        write_gml(&graph, path.to_str().unwrap())?;
        let output = std::fs::read_to_string(path)?;
        assert!(output.contains("label \"0\""));
        assert!(output.contains("centroid \"_networkx_list_start\"\n    centroid \"0_0_1\""));
        assert!(output.contains("members \"_networkx_list_start\"\n    members 0"));
        assert!(output.contains("seqIDs \"_networkx_list_start\"\n    seqIDs \"0_0_1\""));
        assert!(output.contains("dna \"_networkx_list_start\"\n    dna \"ATGAAA\""));
        assert!(output.contains("annotation \"a&#34;b&#38;c\""));
        Ok(())
    }

    #[test]
    fn presence_absence_preserves_genome_order_and_original_ids() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let mut graph = Graph::default();
        graph.add_node(GeneNode {
            annotation: "gene".into(),
            size: 2,
            seq_ids: HashSet::from(["2_0_1".into(), "0_0_2".into(), "0_0_1".into()]),
            ..Default::default()
        });
        let names = HashMap::from([
            (2, "last".into()),
            (0, "first".into()),
            (1, "middle".into()),
        ]);
        let ids = HashMap::from([
            ("0_0_1".into(), "original_a".into()),
            ("0_0_2".into(), "original_b".into()),
            ("2_0_1".into(), "original_c".into()),
        ]);
        generate_roary_gene_presence_absence(&graph, &names, &ids, dir.path().to_str().unwrap())?;
        let csv = std::fs::read_to_string(dir.path().join("gene_presence_absence.csv"))?;
        let lines: Vec<_> = csv.lines().collect();
        assert_eq!(
            lines[0],
            "Gene,Non-unique Gene name,Annotation,first,middle,last"
        );
        assert_eq!(lines[1], "gene,gene,,original_a;original_b,,original_c");
        Ok(())
    }

    #[test]
    fn sparse_rows_preserve_all_output_bytes_and_stable_ties() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let mut graph = Graph::default();
        graph.add_node(GeneNode {
            annotation: "gene".into(),
            description: "first".into(),
            size: 1,
            seq_ids: HashSet::from(["0_0_1".into()]),
            lengths: vec![3],
            ..Default::default()
        });
        graph.add_node(GeneNode {
            annotation: "gene".into(),
            description: "second".into(),
            size: 1,
            seq_ids: HashSet::from(["1_0_1".into(), "1_0_2".into()]),
            lengths: vec![3, 6],
            ..Default::default()
        });
        graph.add_node(GeneNode {
            annotation: "other".into(),
            size: 2,
            seq_ids: HashSet::from(["0_0_2".into(), "1_0_3".into()]),
            lengths: vec![6],
            ..Default::default()
        });
        let names = HashMap::from([(0, "A".into()), (1, "B".into())]);
        let ids = HashMap::from([("1_0_1".into(), "".into()), ("1_0_2".into(), "real".into())]);
        generate_roary_gene_presence_absence(&graph, &names, &ids, dir.path().to_str().unwrap())?;
        assert_eq!(std::fs::read_to_string(dir.path().join("gene_presence_absence.csv"))?,
            "Gene,Non-unique Gene name,Annotation,A,B\nother,other,,0_0_2,1_0_3\ngene,gene,first,0_0_1,\ngroup_0,gene,second,,real\n");
        assert_eq!(
            std::fs::read_to_string(dir.path().join("gene_presence_absence.Rtab"))?,
            "Gene\tA\tB\nother\t1\t1\ngene\t1\t0\ngroup_0\t0\t1\n"
        );
        assert_eq!(std::fs::read_to_string(dir.path().join("gene_presence_absence_roary.csv"))?,
            "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome Fragment,Order within Fragment,Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,Avg group size nuc,A,B\nother,other,,2,2,1.00,3,1,,,,6,6,6.0,0_0_2,1_0_3\ngene,gene,first,1,1,1.00,1,1,,,,3,3,3.0,0_0_1,\ngroup_0,gene,second,1,2,2.00,2,1,,,,3,6,4.5,,real\n");
        Ok(())
    }
}
