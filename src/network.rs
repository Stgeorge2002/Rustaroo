//! Pangenome network generation from CD-HIT clusters
//!

use anyhow::{Context, Result};
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableGraph;
use petgraph::Undirected;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::atomic::{AtomicU64, Ordering};

pub type Graph<N, E> = StableGraph<N, E, Undirected>;

/// Gene node in the pangenome graph
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GeneNode {
    pub size: usize,
    #[serde(skip)]
    pub creation_order: usize,
    pub centroid: Vec<String>,
    pub max_len_id: usize,
    pub members: BTreeSet<usize>,
    pub seq_ids: HashSet<String>,
    pub has_end: bool,
    pub protein: Vec<String>,
    pub dna: Vec<String>,
    pub annotation: String,
    pub description: String,
    pub lengths: Vec<usize>,
    pub long_centroid_id: (usize, String),
    pub paralog: bool,
    pub merged_dna: bool,
}

impl Default for GeneNode {
    fn default() -> Self {
        Self {
            size: 0,
            creation_order: 0,
            centroid: Vec::new(),
            max_len_id: 0,
            members: BTreeSet::new(),
            seq_ids: HashSet::new(),
            has_end: false,
            protein: Vec::new(),
            dna: Vec::new(),
            annotation: String::new(),
            description: String::new(),
            lengths: Vec::new(),
            long_centroid_id: (0, String::new()),
            paralog: false,
            merged_dna: false,
        }
    }
}

/// Edge data in the pangenome graph
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct GeneEdge {
    pub size: usize,
    pub members: BTreeSet<usize>,
    #[serde(skip)]
    pub insertion_order: u64,
}

pub(crate) fn add_ordered_edge(
    graph: &mut Graph<GeneNode, GeneEdge>,
    first: NodeIndex,
    second: NodeIndex,
    mut edge: GeneEdge,
) -> EdgeIndex {
    // Petgraph can reuse edge slots, so a slot is not an insertion timestamp.
    static NEXT_EDGE: AtomicU64 = AtomicU64::new(1);
    edge.insertion_order = NEXT_EDGE.fetch_add(1, Ordering::Relaxed);
    graph.add_edge(first, second, edge)
}

pub(crate) fn ordered_edges(
    graph: &Graph<GeneNode, GeneEdge>,
    node: NodeIndex,
    scratch: &mut Vec<(u64, EdgeIndex, NodeIndex)>,
) {
    use petgraph::visit::EdgeRef;
    scratch.clear();
    scratch.extend(
        graph
            .edges(node)
            .map(|edge| (edge.weight().insertion_order, edge.id(), edge.target())),
    );
    scratch.sort_unstable_by_key(|&(order, edge, _)| (order, edge));
}

#[derive(Debug)]
struct ClusterData {
    prot_sequence: String,
    dna_sequence: String,
    annotation: String,
    description: String,
}

/// Generate pangenome network from CD-HIT output
pub fn generate_network(
    cluster_file: &str,
    data_file: &str,
    prot_seq_file: &str,
    all_dna: bool,
) -> Result<(
    Graph<GeneNode, GeneEdge>,
    HashMap<String, Vec<(NodeIndex, usize)>>,
    HashMap<String, String>,
)> {
    // Parse cluster file
    let (seq_to_cluster, cluster_centroids, cluster_members) = parse_cluster_file(cluster_file)?;

    // Determine paralogs
    let mut paralogs = HashSet::new();
    for (cluster, members) in &cluster_members {
        let genomes: Vec<&str> = members
            .iter()
            .map(|s| s.split('_').next().unwrap())
            .collect();
        let unique_genomes: HashSet<_> = genomes.into_iter().collect();
        if members.len() != unique_genomes.len() {
            paralogs.insert(*cluster);
        }
    }

    // Load sequence metadata
    let cluster_data = load_cluster_data(data_file, &cluster_centroids, &seq_to_cluster)?;

    // Load sequence IDs from protein file
    let seq_ids = load_sequence_ids(prot_seq_file)?;

    // Build graph
    let (graph, centroid_context, seqid_to_centroid) = build_graph(
        &seq_ids,
        &seq_to_cluster,
        &cluster_centroids,
        &cluster_data,
        &paralogs,
        all_dna,
    )?;

    Ok((graph, centroid_context, seqid_to_centroid))
}

/// Parse CD-HIT .clstr file
fn parse_cluster_file(
    cluster_file: &str,
) -> Result<(
    HashMap<String, usize>,
    HashMap<usize, String>,
    HashMap<usize, Vec<String>>,
)> {
    let file = File::open(cluster_file)
        .with_context(|| format!("Failed to open cluster file: {}", cluster_file))?;
    let reader = BufReader::new(file);

    let mut seq_to_cluster = HashMap::new();
    let mut cluster_centroids = HashMap::new();
    let mut cluster_members: HashMap<usize, Vec<String>> = HashMap::new();
    let mut current_cluster = 0;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // New cluster
            current_cluster = line
                .split_whitespace()
                .last()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0);
        } else if line.contains('>') {
            // Sequence in cluster
            let parts: Vec<&str> = line.split('>').collect();
            if parts.len() >= 2 {
                let seq_id = parts[1].split("...").next().unwrap_or("").to_string();
                seq_to_cluster.insert(seq_id.clone(), current_cluster);

                // Get genome ID (first part before _)
                let genome_id = seq_id.split('_').next().unwrap_or("").to_string();
                cluster_members
                    .entry(current_cluster)
                    .or_insert_with(Vec::new)
                    .push(genome_id);

                // Check if this is the representative (marked with *)
                if line.trim().ends_with('*') {
                    cluster_centroids.insert(current_cluster, seq_id);
                }
            }
        }
    }

    Ok((seq_to_cluster, cluster_centroids, cluster_members))
}

/// Load cluster metadata from CSV
fn load_cluster_data(
    data_file: &str,
    cluster_centroids: &HashMap<usize, String>,
    seq_to_cluster: &HashMap<String, usize>,
) -> Result<HashMap<usize, ClusterData>> {
    let file = File::open(data_file)
        .with_context(|| format!("Failed to open data file: {}", data_file))?;
    let reader = BufReader::new(file);

    let centroid_ids: HashSet<_> = cluster_centroids.values().cloned().collect();
    let mut cluster_data = HashMap::new();
    let mut standard_schema = false;

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 {
            standard_schema = line.starts_with("gff_file,scaffold_name,");
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() >= 7 {
            let seq_id = fields[2].to_string();
            if centroid_ids.contains(&seq_id) {
                if let Some(&cluster) = seq_to_cluster.get(&seq_id) {
                    cluster_data.insert(
                        cluster,
                        ClusterData {
                            prot_sequence: fields[4].to_string(),
                            dna_sequence: fields[if standard_schema { 5 } else { 3 }].to_string(),
                            annotation: fields[6].to_string(),
                            description: if standard_schema {
                                fields.get(7).unwrap_or(&"").to_string()
                            } else {
                                fields[5].to_string()
                            },
                        },
                    );
                }
            }
        }
    }

    if cluster_data.len() != cluster_centroids.len() {
        anyhow::bail!(
            "Missing gene metadata for {} cluster centroids",
            cluster_centroids.len() - cluster_data.len()
        );
    }
    Ok(cluster_data)
}

/// Load sequence IDs from FASTA
fn load_sequence_ids(fasta_file: &str) -> Result<Vec<String>> {
    let file = File::open(fasta_file)?;
    let reader = BufReader::new(file);
    let mut seq_ids = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            let id = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            seq_ids.push(id);
        }
    }

    Ok(seq_ids)
}

/// Build the pangenome graph
fn build_graph(
    seq_ids: &[String],
    seq_to_cluster: &HashMap<String, usize>,
    cluster_centroids: &HashMap<usize, String>,
    cluster_data: &HashMap<usize, ClusterData>,
    paralogs: &HashSet<usize>,
    all_dna: bool,
) -> Result<(
    Graph<GeneNode, GeneEdge>,
    HashMap<String, Vec<(NodeIndex, usize)>>,
    HashMap<String, String>,
)> {
    let mut graph: Graph<GeneNode, GeneEdge> = Graph::default();
    let mut centroid_context: HashMap<String, Vec<(NodeIndex, usize)>> = HashMap::new();
    let mut seqid_to_centroid = HashMap::new();
    let mut node_map: HashMap<usize, NodeIndex> = HashMap::new();
    let mut edge_index = HashMap::new();
    let mut prev: Option<NodeIndex> = None;
    let mut previous_location: Option<(usize, usize)> = None;
    // Panaroo keeps the CD-HIT cluster ID for regular nodes and numbers
    // paralog occurrences after the cluster IDs, in encounter order.
    let mut next_paralog_id = cluster_centroids.keys().max().map_or(0, |id| id + 2);

    for seq_id in seq_ids {
        let current_cluster = match seq_to_cluster.get(seq_id) {
            Some(&c) => c,
            None => continue,
        };

        let centroid = match cluster_centroids.get(&current_cluster) {
            Some(c) => c.clone(),
            None => continue,
        };

        seqid_to_centroid.insert(seq_id.clone(), centroid.clone());

        let parts: Vec<&str> = seq_id.split('_').collect();
        let genome_id: usize = parts[0].parse().unwrap_or(0);
        let scaffold = parts
            .get(1)
            .context("Sequence ID is missing its scaffold index")?
            .parse::<usize>()?;
        let is_contig_start = previous_location != Some((genome_id, scaffold));
        previous_location = Some((genome_id, scaffold));

        if is_contig_start {
            // Start of a contig
            if let Some(p) = prev {
                graph[p].has_end = true;
            }

            if node_map.contains_key(&current_cluster) && !paralogs.contains(&current_cluster) {
                // Existing non-paralog node
                let node_idx = node_map[&current_cluster];
                graph[node_idx].size += 1;
                graph[node_idx].members.insert(genome_id);
                graph[node_idx].seq_ids.insert(seq_id.clone());
                graph[node_idx].has_end = true;
                if let Some(data) = cluster_data.get(&current_cluster) {
                    graph[node_idx].lengths.push(data.dna_sequence.len());
                    if all_dna {
                        graph[node_idx].dna.push(data.dna_sequence.clone());
                    }
                }
                prev = Some(node_idx);
            } else {
                // New node (or paralog)
                if paralogs.contains(&current_cluster) {}
                let node_idx = add_node(
                    &mut graph,
                    &mut node_map,
                    current_cluster,
                    cluster_centroids,
                    cluster_data,
                    genome_id,
                    seq_id,
                    paralogs.contains(&current_cluster),
                    true,
                )?;

                if paralogs.contains(&current_cluster) {
                    graph[node_idx].creation_order = next_paralog_id;
                    next_paralog_id += 1;
                    centroid_context
                        .entry(centroid.clone())
                        .or_insert_with(Vec::new)
                        .push((node_idx, genome_id));
                }
                prev = Some(node_idx);
            }
        } else {
            // Middle of contig
            let is_paralog = paralogs.contains(&current_cluster);

            if is_paralog {
                // Create new paralog node
                let node_idx = add_node(
                    &mut graph,
                    &mut node_map,
                    current_cluster,
                    cluster_centroids,
                    cluster_data,
                    genome_id,
                    seq_id,
                    true,
                    false,
                )?;

                graph[node_idx].creation_order = next_paralog_id;
                next_paralog_id += 1;
                centroid_context
                    .entry(centroid.clone())
                    .or_insert_with(Vec::new)
                    .push((node_idx, genome_id));

                if let Some(p) = prev {
                    let edge = GeneEdge {
                        size: 1,
                        members: {
                            let mut s = BTreeSet::new();
                            s.insert(genome_id);
                            s
                        },
                        ..GeneEdge::default()
                    };
                    add_indexed_edge(&mut graph, &mut edge_index, p, node_idx, edge);
                }
                prev = Some(node_idx);
            } else {
                // Non-paralog
                if !node_map.contains_key(&current_cluster) {
                    let node_idx = add_node(
                        &mut graph,
                        &mut node_map,
                        current_cluster,
                        cluster_centroids,
                        cluster_data,
                        genome_id,
                        seq_id,
                        false,
                        false,
                    )?;

                    if let Some(p) = prev {
                        let edge = GeneEdge {
                            size: 1,
                            members: {
                                let mut s = BTreeSet::new();
                                s.insert(genome_id);
                                s
                            },
                            ..GeneEdge::default()
                        };
                        add_indexed_edge(&mut graph, &mut edge_index, p, node_idx, edge);
                    }
                    prev = Some(node_idx);
                } else {
                    let node_idx = node_map[&current_cluster];
                    graph[node_idx].size += 1;
                    graph[node_idx].members.insert(genome_id);
                    graph[node_idx].seq_ids.insert(seq_id.clone());

                    if let Some(data) = cluster_data.get(&current_cluster) {
                        graph[node_idx].lengths.push(data.dna_sequence.len());
                        if all_dna {
                            graph[node_idx].dna.push(data.dna_sequence.clone());
                        }
                    }

                    if let Some(p) = prev {
                        if let Some(edge_idx) = indexed_edge(&edge_index, p, node_idx) {
                            graph[edge_idx].size += 1;
                            graph[edge_idx].members.insert(genome_id);
                        } else {
                            let edge = GeneEdge {
                                size: 1,
                                members: {
                                    let mut s = BTreeSet::new();
                                    s.insert(genome_id);
                                    s
                                },
                                ..GeneEdge::default()
                            };
                            add_indexed_edge(&mut graph, &mut edge_index, p, node_idx, edge);
                        }
                    }
                    prev = Some(node_idx);
                }
            }
        }
    }

    if let Some(p) = prev {
        graph[p].has_end = true;
    }

    Ok((graph, centroid_context, seqid_to_centroid))
}

fn indexed_edge(
    index: &HashMap<(NodeIndex, NodeIndex), EdgeIndex>,
    first: NodeIndex,
    second: NodeIndex,
) -> Option<EdgeIndex> {
    index
        .get(&(first, second))
        .or_else(|| index.get(&(second, first)))
        .copied()
}

fn add_indexed_edge(
    graph: &mut Graph<GeneNode, GeneEdge>,
    index: &mut HashMap<(NodeIndex, NodeIndex), EdgeIndex>,
    first: NodeIndex,
    second: NodeIndex,
    edge: GeneEdge,
) {
    let edge = add_ordered_edge(graph, first, second, edge);
    index.insert((first, second), edge);
}

fn add_node(
    graph: &mut Graph<GeneNode, GeneEdge>,
    node_map: &mut HashMap<usize, NodeIndex>,
    cluster_id: usize,
    cluster_centroids: &HashMap<usize, String>,
    cluster_data: &HashMap<usize, ClusterData>,
    genome_id: usize,
    seq_id: &str,
    is_paralog: bool,
    has_end: bool,
) -> Result<NodeIndex> {
    let centroid = cluster_centroids
        .get(&cluster_id)
        .map(|s| s.clone())
        .unwrap_or_default();
    let data = cluster_data.get(&cluster_id);

    let node = GeneNode {
        size: 1,
        creation_order: cluster_id,
        centroid: vec![centroid.clone()],
        max_len_id: 0,
        members: {
            let mut s = BTreeSet::new();
            s.insert(genome_id);
            s
        },
        seq_ids: {
            let mut s = HashSet::new();
            s.insert(seq_id.to_string());
            s
        },
        has_end,
        protein: data
            .as_ref()
            .map(|d| vec![d.prot_sequence.clone()])
            .unwrap_or_default(),
        dna: data
            .as_ref()
            .map(|d| vec![d.dna_sequence.clone()])
            .unwrap_or_default(),
        annotation: data
            .as_ref()
            .map(|d| d.annotation.clone())
            .unwrap_or_default(),
        description: data
            .as_ref()
            .map(|d| d.description.clone())
            .unwrap_or_default(),
        lengths: data
            .as_ref()
            .map(|d| vec![d.dna_sequence.len()])
            .unwrap_or_default(),
        long_centroid_id: (
            data.as_ref().map(|d| d.dna_sequence.len()).unwrap_or(0),
            centroid,
        ),
        paralog: is_paralog,
        merged_dna: false,
    };

    let node_idx = graph.add_node(node);
    if !is_paralog {
        node_map.insert(cluster_id, node_idx);
    }
    Ok(node_idx)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexed_edges_match_graph_lookup_with_parallel_edges_and_loops() {
        let mut graph = Graph::default();
        let a = graph.add_node(GeneNode::default());
        let b = graph.add_node(GeneNode::default());
        let c = graph.add_node(GeneNode::default());
        let mut index = HashMap::new();
        for (first, second) in [(a, b), (b, a), (a, a), (b, c), (a, b), (a, a), (b, a)] {
            add_indexed_edge(
                &mut graph,
                &mut index,
                first,
                second,
                GeneEdge {
                    size: 1,
                    members: BTreeSet::new(),
                    ..GeneEdge::default()
                },
            );
            for x in [a, b, c] {
                for y in [a, b, c] {
                    assert_eq!(indexed_edge(&index, x, y), graph.find_edge(x, y));
                }
            }
        }
    }
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_gene_node_creation() {
        let node = GeneNode::default();
        assert_eq!(node.size, 0);
        assert!(!node.paralog);
    }

    #[test]
    fn test_load_cluster_data_from_seven_column_csv() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "clustering_id,gff_file,annotation_id,dna_sequence,prot_sequence,description,name"
        )
        .unwrap();
        writeln!(file, "0,genome,0_1,ATG,MP,description,gene").unwrap();

        let mut centroids = HashMap::new();
        centroids.insert(0, "0_1".to_string());
        let mut seq_to_cluster = HashMap::new();
        seq_to_cluster.insert("0_1".to_string(), 0);

        let data =
            load_cluster_data(file.path().to_str().unwrap(), &centroids, &seq_to_cluster).unwrap();

        let cluster = data.get(&0).unwrap();
        assert_eq!(cluster.dna_sequence, "ATG");
        assert_eq!(cluster.prot_sequence, "MP");
        assert_eq!(cluster.description, "description");
        assert_eq!(cluster.annotation, "gene");
    }
}
