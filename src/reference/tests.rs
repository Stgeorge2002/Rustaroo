use super::*;
use crate::output::{generate_roary_gene_presence_absence, write_gml};
use genetic_code::back_translate;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use std::collections::BTreeSet;

fn protein(length: usize) -> Vec<u8> {
    let mut state = 0x1234_5678_9abc_def0u64;
    let mut protein: Vec<_> = (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACDEFGHIKLMNPQRSTVWY"[state as usize % 20]
        })
        .collect();
    protein[0] = b'M';
    protein
}

fn snapshot(graph: &Graph<GeneNode, GeneEdge>) -> Result<String> {
    Ok(serde_json::to_string(&(
        graph
            .node_indices()
            .map(|index| (index.index(), &graph[index]))
            .collect::<Vec<_>>(),
        graph
            .edge_references()
            .map(|edge| (edge.source().index(), edge.target().index(), edge.weight()))
            .collect::<Vec<_>>(),
    ))?)
}

fn metadata_header() -> &'static str {
    "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description\n"
}

#[test]
fn all_members_not_just_centroids_choose_references_without_changing_pangenome() -> Result<()> {
    let dir = tempfile::TempDir::new()?;
    let full = back_translate(&protein(100));
    let short = back_translate(&protein(100)[..70]);
    let other = "ATGAAATAG";
    let mut metadata = metadata_header().to_string();
    let mut inputs = Vec::new();
    for (sample, dna) in [&short, &full].iter().enumerate() {
        let start = if sample == 0 { 1 } else { 201 };
        let end = start + dna.len() - 1;
        let gff = dir.path().join(format!("sample{sample}.gff"));
        let fasta = format!("{}{}{}", "C".repeat(start - 1), dna, "C".repeat(300));
        std::fs::write(&gff, format!(
            "##gff-version 3\ncontig\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tID=a\n##FASTA\n>contig\n{fasta}\n",
        ))?;
        inputs.push(gff.to_str().unwrap().to_string());
        metadata.push_str(&format!(
            "sample{sample},contig,{sample}_0_0,a,,{dna},same,\n"
        ));
    }
    // This family shares the other family's centroid, but is a distinct final
    // paralog node. Its own accepted allele must still produce a reference.
    metadata.push_str(&format!(
        "sample0,contig,0_refound_0,0_refound_0,,{other},refound,location:0:9;strand:+\n"
    ));
    let data = dir.path().join("gene_data.csv");
    std::fs::write(&data, metadata)?;
    let mut graph = Graph::default();
    let first = graph.add_node(GeneNode {
        annotation: "same".into(),
        centroid: vec!["0_0_0".into()],
        dna: vec![short.clone()],
        protein: vec!["centroid_only".into()],
        seq_ids: HashSet::from(["0_0_0".into(), "1_0_0".into()]),
        members: BTreeSet::from([0, 1]),
        size: 2,
        ..Default::default()
    });
    let second = graph.add_node(GeneNode {
        annotation: "same".into(),
        centroid: vec!["0_0_0".into()],
        seq_ids: HashSet::from(["0_refound_0".into()]),
        members: BTreeSet::from([0]),
        size: 1,
        ..Default::default()
    });
    graph.add_edge(
        first,
        second,
        GeneEdge {
            size: 1,
            members: BTreeSet::from([0]),
            ..GeneEdge::default()
        },
    );
    let before = snapshot(&graph)?;
    let names = generate_roary_gene_presence_absence(
        &graph,
        &HashMap::from([(0, "sample0".into()), (1, "sample1".into())]),
        &HashMap::new(),
        dir.path().to_str().unwrap(),
    )?;
    let matrix_paths = [
        "gene_presence_absence.csv",
        "gene_presence_absence_roary.csv",
        "gene_presence_absence.Rtab",
    ];
    let matrices: Vec<_> = matrix_paths
        .iter()
        .map(|name| std::fs::read(dir.path().join(name)))
        .collect::<std::io::Result<_>>()?;
    generate_pan_genome_reference(
        &graph,
        data.to_str().unwrap(),
        &inputs,
        &names,
        dir.path().to_str().unwrap(),
    )?;
    let fasta = std::fs::read_to_string(dir.path().join("pan_genome_reference.fa"))?;
    let expected_full = full
        .as_bytes()
        .chunks(60)
        .map(|chunk| std::str::from_utf8(chunk).unwrap())
        .collect::<Vec<_>>()
        .join("\n");
    assert_eq!(
        fasta,
        format!(">same\n{expected_full}\n>group_0\n{other}\n")
    );
    let report = std::fs::read(dir.path().join("reference_report.tsv"))?;
    for threads in [1, 4] {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()?
            .install(|| {
                generate_pan_genome_reference(
                    &graph,
                    data.to_str().unwrap(),
                    &inputs,
                    &names,
                    dir.path().to_str().unwrap(),
                )
            })?;
        assert_eq!(
            std::fs::read_to_string(dir.path().join("pan_genome_reference.fa"))?,
            fasta
        );
        assert_eq!(
            std::fs::read(dir.path().join("reference_report.tsv"))?,
            report
        );
    }
    assert_eq!(before, snapshot(&graph)?);
    for (path, expected) in matrix_paths.iter().zip(matrices) {
        assert_eq!(std::fs::read(dir.path().join(path))?, expected);
    }
    Ok(())
}

fn call(pools: &mut Pools, dna: String, sample: u32, start: u32, strand: i8) -> Call {
    let protein = translate(&dna);
    let dna_id = pools
        .dna
        .iter()
        .position(|existing| existing == &dna)
        .unwrap_or(pools.dna.len());
    if dna_id == pools.dna.len() {
        pools.dna.push(dna);
    }
    let protein_id = pools
        .protein
        .iter()
        .position(|existing| existing == &protein)
        .unwrap_or(pools.protein.len());
    if protein_id == pools.protein.len() {
        pools.protein.push(protein);
    }
    Call {
        sample,
        contig: 0,
        start,
        end: start + pools.dna[dna_id].len() as u32 - 1,
        strand,
        dna: dna_id as u32,
        protein: protein_id as u32,
        cluster: 0,
        locus: None,
        refound: false,
    }
}

#[test]
fn refound_short_copies_stay_members_but_cannot_replace_annotated_reference() -> Result<()> {
    let full = protein(100);
    let mut pools = Pools::default();
    let mut calls = Vec::new();
    for sample in 0..12 {
        let mut member = call(
            &mut pools,
            back_translate(if sample < 2 { &full } else { &full[..80] }),
            sample,
            201,
            1,
        );
        member.refound = sample >= 2;
        calls.push(member);
    }
    let indexes: Vec<_> = (0..calls.len()).collect();
    let loci: Vec<_> = (0..calls.len() as u32).collect();
    let lengths = vec![vec![4000]; 12];
    let alleles = alleles_of(&indexes, &calls, &pools, &lengths, &loci);
    let short = alleles
        .iter()
        .find(|allele| allele.dna == back_translate(&full[..80]))
        .unwrap();
    assert_eq!((short.count, short.finished), (10, 0));
    let choice = select(
        0,
        &indexes,
        &calls,
        &pools,
        &OtherGenes::new(FxMap::default()),
        &loci,
        &lengths,
        &(0..calls.len())
            .map(|i| format!("{i}_0_0"))
            .collect::<Vec<_>>(),
    )?;
    assert_eq!(choice.choice.dna, back_translate(&full));
    assert_eq!(choice.exact_samples, 2);
    assert_eq!(choice.excluded_refound_short, 1);
    assert_eq!(calls.len(), 12);
    Ok(())
}

#[test]
fn distinct_sample_support_and_symmetric_contig_edges_match_rustaroo() {
    let mut pools = Pools::default();
    let dna = back_translate(&protein(30));
    let calls = vec![
        call(&mut pools, dna.clone(), 0, 201, 1),
        call(&mut pools, dna.clone(), 0, 401, 1),
        call(&mut pools, dna.clone(), 1, 91, 1),
        call(&mut pools, dna, 2, 201, -1),
    ];
    let lengths = vec![vec![4000], vec![4000], vec![calls[3].end as usize + 90]];
    let alleles = alleles_of(&[0, 1, 2, 3], &calls, &pools, &lengths, &[0, 1, 2, 3]);
    assert_eq!(alleles.len(), 1);
    assert_eq!(alleles[0].count, 3);
    // Both copies exactly 90 bases from an end are clipped; duplicate copies
    // in the internal sample contribute one finished-end vote.
    assert_eq!(alleles[0].finished, 1);
}

#[test]
fn fragment_models_use_one_existing_family_locus_on_either_strand() -> Result<()> {
    let full = protein(80);
    for reverse in [false, true] {
        let mut pools = Pools::default();
        let (head, tail) = if reverse {
            (&full[40..], &full[..40])
        } else {
            (&full[..40], &full[40..])
        };
        let strand = if reverse { -1 } else { 1 };
        let a = call(&mut pools, back_translate(head), 0, 201, strand);
        let b = call(&mut pools, back_translate(tail), 0, a.end + 1, strand);
        let calls = vec![a, b];
        let loci = membership_loci(&calls, &pools);
        assert_eq!(loci[0], loci[1]);
        assert_eq!(
            joined_allele(&[0, 1], &calls, &pools, &loci),
            Some(back_translate(&full))
        );
        assert!(joined_allele(&[0, 1], &calls, &pools, &[0, 1]).is_none());
        let selected = select(
            0,
            &[0, 1],
            &calls,
            &pools,
            &OtherGenes::new(FxMap::default()),
            &loci,
            &[vec![4000]],
            &["head".into(), "tail".into()],
        )?;
        assert_eq!(selected.choice.selection, "joined_fragments");
        assert_eq!(selected.choice.dna, back_translate(&full));
        assert_eq!(selected.exact_samples, 0);
        assert_eq!(selected.excluded_refound_short, 0);
        assert_eq!(selected.sources, ["head", "tail"]);
    }
    for separate_contig in [false, true] {
        let mut pools = Pools::default();
        let a = call(&mut pools, back_translate(&full), 0, 201, 1);
        let mut b = call(&mut pools, back_translate(&full), 0, a.end + 1, 1);
        b.contig = u32::from(separate_contig);
        let calls = vec![a, b];
        let loci = membership_loci(&calls, &pools);
        assert!(joined_allele(&[0, 1], &calls, &pools, &loci).is_none());
    }
    Ok(())
}

#[test]
fn refound_half_open_coordinates_and_missing_evidence_are_explicit() -> Result<()> {
    for separator in [":", "-"] {
        let location =
            refound_location(&format!("location:0{separator}9;strand:-"), "0_refound_0")?;
        assert_eq!((location.start, location.end, location.strand), (1, 9, -1));
    }
    assert!(refound_location("location:0:9;strand:?", "0_refound_0").is_err());
    let dir = tempfile::TempDir::new()?;
    let data = dir.path().join("gene_data.csv");
    std::fs::write(&data, metadata_header())?;
    let mut graph = Graph::default();
    let node = graph.add_node(GeneNode {
        seq_ids: HashSet::from(["0_0_0".into()]),
        ..Default::default()
    });
    let error = generate_pan_genome_reference(
        &graph,
        data.to_str().unwrap(),
        &["unused.gff".into()],
        &HashMap::from([(node, "gene".into())]),
        dir.path().to_str().unwrap(),
    )
    .unwrap_err();
    assert!(error.to_string().contains("Missing metadata"));
    assert!(!dir.path().join("pan_genome_reference.fa").exists());
    Ok(())
}

#[test]
#[ignore = "requires RUSTAROO_REFERENCE_FIXTURE containing nodes, edges, gene_data_file, input_files"]
fn finalized_real_graph_remains_unchanged() -> Result<()> {
    #[derive(serde::Deserialize)]
    struct Fixture {
        nodes: Vec<(usize, GeneNode)>,
        edges: Vec<(usize, usize, GeneEdge)>,
        gene_data_file: String,
        input_files: Vec<String>,
    }
    let fixture: Fixture = serde_json::from_reader(BufReader::new(File::open(std::env::var(
        "RUSTAROO_REFERENCE_FIXTURE",
    )?)?))?;
    let mut graph = Graph::default();
    let mut ids = HashMap::new();
    for (id, node) in fixture.nodes {
        ids.insert(id, graph.add_node(node));
    }
    for (first, second, edge) in fixture.edges {
        graph.add_edge(ids[&first], ids[&second], edge);
    }
    let dir = tempfile::TempDir::new()?;
    let before = snapshot(&graph)?;
    let isolates = fixture
        .input_files
        .iter()
        .enumerate()
        .map(|(i, file)| {
            (
                i,
                Path::new(file)
                    .file_stem()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .to_string(),
            )
        })
        .collect();
    let names = generate_roary_gene_presence_absence(
        &graph,
        &isolates,
        &HashMap::new(),
        dir.path().to_str().unwrap(),
    )?;
    let paths = [
        "gene_presence_absence.csv",
        "gene_presence_absence_roary.csv",
        "gene_presence_absence.Rtab",
    ];
    let original: Vec<_> = paths
        .iter()
        .map(|name| std::fs::read(dir.path().join(name)))
        .collect::<std::io::Result<_>>()?;
    write_gml(&graph, dir.path().join("before.gml").to_str().unwrap())?;
    let started = std::time::Instant::now();
    generate_pan_genome_reference(
        &graph,
        &fixture.gene_data_file,
        &fixture.input_files,
        &names,
        dir.path().to_str().unwrap(),
    )?;
    println!(
        "Reference selection: {:.3}s for {} unchanged families",
        started.elapsed().as_secs_f64(),
        graph.node_count()
    );
    assert_eq!(snapshot(&graph)?, before);
    write_gml(&graph, dir.path().join("after.gml").to_str().unwrap())?;
    assert_eq!(
        std::fs::read(dir.path().join("before.gml"))?,
        std::fs::read(dir.path().join("after.gml"))?
    );
    for (path, original) in paths.iter().zip(original) {
        assert_eq!(std::fs::read(dir.path().join(path))?, original);
    }
    let fasta = std::fs::read_to_string(dir.path().join("pan_genome_reference.fa"))?;
    assert_eq!(
        fasta.lines().filter(|line| line.starts_with('>')).count(),
        graph.node_count()
    );
    Ok(())
}
