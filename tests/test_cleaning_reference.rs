use anyhow::Result;
use petgraph::graph::NodeIndex;
use rustaroo::clean_network::collapse_families_indexed;
use rustaroo::family_index::FamilyIndex;
use rustaroo::merge_nodes::merge_node_cluster;
use rustaroo::network::{GeneEdge, GeneNode, Graph};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::time::Instant;

#[derive(Deserialize)]
struct Fixture {
    nodes: Vec<(usize, GeneNode)>,
    edges: Vec<(usize, usize, GeneEdge)>,
    seqid_to_centroid: HashMap<String, String>,
}

#[derive(Serialize)]
struct Cleaned {
    partitions: BTreeSet<BTreeSet<String>>,
    edges: BTreeSet<(BTreeSet<String>, BTreeSet<String>, BTreeSet<usize>)>,
    seconds: f64,
}

#[test]
fn merging_isolated_nodes_removes_originals_without_moving_live_nodes() -> Result<()> {
    let mut graph = Graph::default();
    let first = graph.add_node(GeneNode::default());
    let second = graph.add_node(GeneNode::default());
    let live = graph.add_node(GeneNode {
        annotation: "live".into(),
        ..GeneNode::default()
    });
    let actual = merge_node_cluster(&mut graph, &[first, second], NodeIndex::new(99), true)?;
    assert!(graph.node_weight(actual).is_some());
    assert!(graph.node_weight(first).is_none());
    assert!(graph.node_weight(second).is_none());
    assert_eq!(graph[live].annotation, "live");
    assert_eq!(graph.node_count(), 2);
    Ok(())
}

#[test]
fn merging_distinct_centroids_preserves_identical_protein_associations() -> Result<()> {
    let mut graph = Graph::default();
    let first = graph.add_node(GeneNode {
        centroid: vec!["a".into()],
        dna: vec!["ATGAAA".into()],
        protein: vec!["MK".into()],
        ..Default::default()
    });
    let second = graph.add_node(GeneNode {
        centroid: vec!["b".into()],
        dna: vec!["ATGAAG".into()],
        protein: vec!["MK".into()],
        ..Default::default()
    });
    let merged = merge_node_cluster(&mut graph, &[first, second], NodeIndex::new(99), false)?;
    assert_eq!(graph[merged].centroid, ["a", "b"]);
    assert_eq!(graph[merged].dna, ["ATGAAA", "ATGAAG"]);
    assert_eq!(graph[merged].protein, ["MK", "MK"]);
    let duplicate = graph.add_node(graph[merged].clone());
    let merged = merge_node_cluster(&mut graph, &[merged, duplicate], NodeIndex::new(100), true)?;
    assert_eq!(graph[merged].centroid, ["a", "b"]);
    assert_eq!(graph[merged].protein, ["MK", "MK"]);
    Ok(())
}

#[test]
fn high_degree_nonmerging_cleaning_terminates_and_preserves_graph() -> Result<()> {
    let mut graph = Graph::default();
    let mut names = HashMap::new();
    let mut mapping = HashMap::new();
    for i in 0..5 {
        let name = format!("{}_0_0", i);
        names.insert(name.clone(), i);
        mapping.insert(name.clone(), name.clone());
        graph.add_node(GeneNode {
            centroid: vec![name.clone()],
            seq_ids: HashSet::from([name]),
            members: BTreeSet::from([i]),
            size: 1,
            ..GeneNode::default()
        });
    }
    for i in 1..5 {
        graph.add_edge(
            NodeIndex::new(0),
            NodeIndex::new(i),
            GeneEdge {
                size: 0,
                members: BTreeSet::new(),
                ..GeneEdge::default()
            },
        );
    }
    let index = FamilyIndex {
        centroid_to_index: names,
        close: vec![Default::default(); 5],
    };
    collapse_families_indexed(&mut graph, &mapping, &index, &[1, 2, 3], true, false, 0.1)?;
    assert_eq!(graph.node_count(), 5);
    assert_eq!(graph.edge_count(), 4);
    Ok(())
}

#[test]
#[ignore = "requires RUSTAROO_REFERENCE_INPUT and RUSTAROO_REFERENCE_OUTPUT fixture paths"]
fn benchmark_reference_graph() -> Result<()> {
    let fixture: Fixture = serde_json::from_reader(std::io::BufReader::new(std::fs::File::open(
        std::env::var("RUSTAROO_REFERENCE_INPUT")?,
    )?))?;
    let mut graph = Graph::default();
    let mut indices = HashMap::new();
    for (id, node) in fixture.nodes {
        indices.insert(id, graph.add_node(node));
    }
    for (first, second, edge) in fixture.edges {
        graph.add_edge(indices[&first], indices[&second], edge);
    }
    let outdir = tempfile::TempDir::new()?;
    let threads = std::env::var("RUSTAROO_REFERENCE_THREADS")
        .unwrap_or_else(|_| "20".into())
        .parse()?;
    let started = Instant::now();
    println!("Building DNA index...");
    let dna = FamilyIndex::build(
        &graph,
        outdir.path().to_str().unwrap(),
        0.98,
        0.0,
        threads,
        true,
    )?;
    println!("Cleaning DNA...");
    collapse_families_indexed(
        &mut graph,
        &fixture.seqid_to_centroid,
        &dna,
        &[1, 2, 3],
        true,
        true,
        0.1,
    )?;
    println!("Building protein index...");
    let protein = FamilyIndex::build(
        &graph,
        outdir.path().to_str().unwrap(),
        0.7,
        0.0,
        threads,
        false,
    )?;
    collapse_families_indexed(
        &mut graph,
        &fixture.seqid_to_centroid,
        &protein,
        &[1, 2, 3],
        true,
        false,
        0.1,
    )?;
    let seconds = started.elapsed().as_secs_f64();
    let partitions = graph
        .node_weights()
        .map(|node| node.seq_ids.iter().cloned().collect())
        .collect();
    let edges = graph
        .edge_indices()
        .map(|edge| {
            let (first, second) = graph.edge_endpoints(edge).unwrap();
            let (mut first, mut second): (BTreeSet<_>, BTreeSet<_>) = (
                graph[first].seq_ids.iter().cloned().collect(),
                graph[second].seq_ids.iter().cloned().collect(),
            );
            if second < first {
                std::mem::swap(&mut first, &mut second);
            }
            (first, second, graph[edge].members.clone())
        })
        .collect();
    let mut writer = std::io::BufWriter::new(std::fs::File::create(std::env::var(
        "RUSTAROO_REFERENCE_OUTPUT",
    )?)?);
    serde_json::to_writer(
        &mut writer,
        &Cleaned {
            partitions,
            edges,
            seconds,
        },
    )?;
    std::io::Write::flush(&mut writer)?;
    println!("Reference graph cleaning: {seconds:.3}s");
    Ok(())
}
