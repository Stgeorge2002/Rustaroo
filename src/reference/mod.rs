//! Rustaroo reference selection over immutable, finalized Rustaroo families.
//!
//! The sequence policy and evidence helpers are ported from Rustaroo revision
//! 6b09ae5a5a401003b6a24195132c15524109e965. See LICENSE.rustaroo.
//! Nothing in this module constructs, cleans, filters, or mutates a pangenome.

mod align;
mod evidence;
pub(crate) mod fxhash;
mod genetic_code;
mod model;
mod reference;
mod synteny;
#[cfg(test)]
mod tests;

use self::evidence::{alleles_of, cluster_reps, joined_allele, joined_models, membership_loci};
use self::fxhash::{HashMap as FxMap, HashSet as FxSet};
use self::genetic_code::{is_intact, translate, uppercase_dna};
use self::model::{Call, Pools};
use self::reference::{choose_reference, Allele, Choice, OtherGenes};
use crate::line_reader::LineReader;
use crate::network::{GeneEdge, GeneNode, Graph};
use anyhow::{bail, Context, Result};
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
#[cfg(test)]
use std::io::BufReader;
use std::path::Path;

pub struct Member {
    id: String,
    sample: usize,
    cluster: u32,
    scaffold: String,
    annotation: String,
    dna: String,
    refound: bool,
    location: Option<Location>,
}

#[derive(Clone, Copy)]
struct Location {
    start: u32,
    end: u32,
    strand: i8,
}

// The maps below are only ever queried, never iterated, so their hasher
// cannot affect any output (see `fxhash`).
struct SampleLocations {
    contigs: FxMap<String, u32>,
    lengths: Vec<usize>,
    annotations: FxMap<(String, String), Location>,
}

/// Annotations needed from one sample: scaffold -> annotation IDs.
type NeededAnnotations = FxMap<String, FxSet<String>>;

struct Selected {
    choice: Choice,
    excluded_refound_short: usize,
    distinct: usize,
    intact: usize,
    exact_samples: usize,
    sources: Vec<String>,
}

/// Replace centroid-only references, without modifying the graph or its members.
pub fn generate_pan_genome_reference(
    graph: &Graph<GeneNode, GeneEdge>,
    gene_data_file: &str,
    input_files: &[String],
    family_names: &HashMap<NodeIndex, String>,
    output_dir: &str,
) -> Result<()> {
    let preload = preload_reference(graph, gene_data_file, input_files.len());
    generate_pan_genome_reference_with(preload, graph, input_files, family_names, output_dir)
}

/// The first family-structure problem found while indexing members, in the
/// original loop order (formatted later, as some messages need the name).
enum StructureError {
    TooManyFamilies(std::num::TryFromIntError),
    NoMembers,
    TooManyMembers(std::num::TryFromIntError),
    Duplicate(String),
}

/// The part of reference generation that does not need family names:
/// the member index of every final family and the member metadata read
/// from gene_data.csv. It can run while the presence/absence output (which
/// produces the names) is generated; `generate_pan_genome_reference_with`
/// then reports any error in exactly the original order.
pub struct ReferencePreload<'g> {
    nodes: Vec<NodeIndex>,
    /// Node position and kind of the first structural error, if any; the
    /// families map and members are only complete without one.
    structure_error: Option<(usize, StructureError)>,
    members: Option<Result<Vec<Member>>>,
    _graph: std::marker::PhantomData<&'g ()>,
}

pub fn preload_reference<'g>(
    graph: &'g Graph<GeneNode, GeneEdge>,
    gene_data_file: &str,
    sample_count: usize,
) -> ReferencePreload<'g> {
    let nodes: Vec<_> = graph.node_indices().collect();
    // Sequence ID -> (family, slot). Slots let `read_members` track which IDs
    // it has seen without a second string set.
    let mut families = FxMap::default();
    let mut structure_error = None;
    'nodes: for (position, &node) in nodes.iter().enumerate() {
        let cluster = match u32::try_from(position) {
            Ok(cluster) => cluster,
            Err(error) => {
                structure_error = Some((position, StructureError::TooManyFamilies(error)));
                break;
            }
        };
        if graph[node].seq_ids.is_empty() {
            structure_error = Some((position, StructureError::NoMembers));
            break;
        }
        for id in &graph[node].seq_ids {
            let slot = match u32::try_from(families.len()) {
                Ok(slot) => slot,
                Err(error) => {
                    structure_error = Some((position, StructureError::TooManyMembers(error)));
                    break 'nodes;
                }
            };
            if families.insert(id.as_str(), (cluster, slot)).is_some() {
                structure_error = Some((position, StructureError::Duplicate(id.clone())));
                break 'nodes;
            }
        }
    }
    let members = structure_error
        .is_none()
        .then(|| read_members(gene_data_file, &families, sample_count));
    ReferencePreload {
        nodes,
        structure_error,
        members,
        _graph: std::marker::PhantomData,
    }
}

/// Finish reference generation from a preload of the same graph.
pub fn generate_pan_genome_reference_with(
    preload: ReferencePreload<'_>,
    graph: &Graph<GeneNode, GeneEdge>,
    input_files: &[String],
    family_names: &HashMap<NodeIndex, String>,
    output_dir: &str,
) -> Result<()> {
    let ReferencePreload {
        nodes,
        structure_error,
        members,
        ..
    } = preload;
    // Same checks, in the same order, as the original single loop: per node,
    // the name checks come before that node's structural checks.
    let mut names = HashSet::new();
    for (position, &node) in nodes.iter().enumerate() {
        let name = family_names
            .get(&node)
            .context("Missing output family name")?;
        if name.is_empty() || name.chars().any(char::is_whitespace) || !names.insert(name) {
            bail!("Reference FASTA requires a unique, nonempty family identifier: {name:?}");
        }
        if let Some((error_position, error)) = &structure_error {
            if *error_position == position {
                match error {
                    StructureError::TooManyFamilies(error) => {
                        return Err(anyhow::Error::new(*error).context("Too many reference families"))
                    }
                    StructureError::NoMembers => {
                        bail!("Family {name} has no member sequences for reference selection")
                    }
                    StructureError::TooManyMembers(error) => {
                        return Err(anyhow::Error::new(*error).context("Too many reference members"))
                    }
                    StructureError::Duplicate(id) => {
                        bail!("Sequence {id} belongs to more than one final family")
                    }
                }
            }
        }
    }

    let mut members = members.expect("members are read when the structure is valid")?;
    let mut needed = vec![NeededAnnotations::default(); input_files.len()];
    for member in &members {
        if !member.refound {
            let scaffolds = &mut needed[member.sample];
            match scaffolds.get_mut(member.scaffold.as_str()) {
                Some(annotations) => {
                    annotations.insert(member.annotation.clone());
                }
                None => {
                    scaffolds.insert(
                        member.scaffold.clone(),
                        FxSet::from_iter([member.annotation.clone()]),
                    );
                }
            }
        }
    }
    let locations: Vec<_> = input_files
        .par_iter()
        .zip(&needed)
        .map(|(path, needed)| read_locations(path, needed))
        .collect::<Result<_>>()?;
    let contig_lengths: Vec<_> = locations
        .iter()
        .map(|sample| sample.lengths.clone())
        .collect();
    let mut pools = Pools::default();
    // Interning maps: IDs are assigned in member order, never by map order.
    let mut dna_ids = FxMap::<String, usize>::default();
    let mut protein_ids = FxMap::<Vec<u8>, usize>::default();
    let mut calls = Vec::with_capacity(members.len());
    // Rustaroo's physical-locus tests require genomic order, not HashSet order.
    // (Stable sorts, as before; the parallel one gives the same order.)
    members.par_sort_by(|a, b| a.id.cmp(&b.id));
    for member in members {
        let sample = &locations[member.sample];
        let contig = *sample.contigs.get(&member.scaffold).with_context(|| {
            format!(
                "Missing contig {} for member {}",
                member.scaffold, member.id
            )
        })?;
        let location = match member.location {
            Some(location) => location,
            None => *sample
                .annotations
                .get(&(member.scaffold, member.annotation))
                .with_context(|| format!("Missing GFF location for member {}", member.id))?,
        };
        if location.start == 0
            || location.end < location.start
            || location.end as usize > sample.lengths[contig as usize]
        {
            bail!("Out-of-range reference evidence for member {}", member.id);
        }
        let dna = match dna_ids.get(member.dna.as_str()) {
            Some(&id) => id,
            None => {
                let id = pools.dna.len();
                pools.dna.push(member.dna.clone());
                dna_ids.insert(member.dna, id);
                id
            }
        };
        calls.push((
            member.id,
            Call {
                sample: u32::try_from(member.sample)?,
                contig,
                start: location.start,
                end: location.end,
                strand: location.strand,
                dna: u32::try_from(dna)?,
                // Filled in below, once every allele is known.
                protein: 0,
                cluster: member.cluster,
                locus: None,
                refound: member.refound,
            },
        ));
    }
    // Protein IDs are assigned in order of each allele's first appearance,
    // which is DNA ID order; translating the alleles in parallel and then
    // interning them in that order gives the same IDs as translating each
    // new allele inside the loop above. (A protein ID never exceeds its DNA
    // ID, so its u32 conversion cannot fail where the DNA one did not.)
    let translated: Vec<Vec<u8>> = pools.dna.par_iter().map(|dna| translate(dna)).collect();
    let mut dna_proteins = Vec::with_capacity(translated.len());
    for protein in translated {
        let id = match protein_ids.get(protein.as_slice()) {
            Some(&id) => id,
            None => {
                let id = pools.protein.len();
                pools.protein.push(protein.clone());
                protein_ids.insert(protein, id);
                id
            }
        };
        dna_proteins.push(u32::try_from(id)?);
    }
    for (_, call) in &mut calls {
        call.protein = dna_proteins[call.dna as usize];
    }
    calls.par_sort_by_key(|(_, call)| (call.sample, call.contig, call.start, call.end));
    let (member_ids, calls): (Vec<_>, Vec<_>) = calls.into_iter().unzip();
    drop(dna_ids);
    drop(protein_ids);
    let mut groups = FxMap::<u32, Vec<usize>>::default();
    for (index, call) in calls.iter().enumerate() {
        groups.entry(call.cluster).or_default().push(index);
    }
    let others = OtherGenes::new(cluster_reps(&groups, &calls, &pools));
    let loci = membership_loci(&calls, &pools);
    let selected: Vec<_> = (0..nodes.len())
        .into_par_iter()
        .map(|cluster| {
            let cluster = u32::try_from(cluster)?;
            select(
                cluster,
                &groups[&cluster],
                &calls,
                &pools,
                &others,
                &loci,
                &contig_lengths,
                &member_ids,
            )
        })
        .collect::<Result<_>>()?;
    let excluded_refound_short: usize = selected
        .iter()
        .map(|selected| selected.excluded_refound_short)
        .sum();

    let mut fasta = BufWriter::new(File::create(
        Path::new(output_dir).join("pan_genome_reference.fa"),
    )?);
    let mut report = BufWriter::new(File::create(
        Path::new(output_dir).join("reference_report.tsv"),
    )?);
    writeln!(report, "gene\treference_bp\tintact\tselection\tsamples_present\tmodal_end_bp\tmodal_end_samples\textension_samples\tsamples_with_exact_reference\tdistinct_alleles\tintact_alleles\tsource_seq_ids")?;
    for (&node, selected) in nodes.iter().zip(selected) {
        let name = &family_names[&node];
        writeln!(fasta, ">{name}")?;
        for chunk in selected.choice.dna.as_bytes().chunks(60) {
            fasta.write_all(chunk)?;
            writeln!(fasta)?;
        }
        let choice = &selected.choice;
        writeln!(
            report,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            name,
            choice.dna.len(),
            if choice.intact { "yes" } else { "no" },
            choice.selection,
            graph[node].members.len(),
            choice.modal_end_bp,
            choice.modal_end_samples,
            choice.extension_samples,
            selected.exact_samples,
            selected.distinct,
            selected.intact,
            selected.sources.join(";"),
        )?;
    }
    fasta.flush()?;
    report.flush()?;
    log::info!(
        "Reference selection: excluded {} distinct short refound candidates",
        excluded_refound_short
    );
    Ok(())
}

fn select(
    cluster: u32,
    indexes: &[usize],
    calls: &[Call],
    pools: &Pools,
    others: &OtherGenes,
    locus_ids: &[u32],
    contig_lengths: &[Vec<usize>],
    member_ids: &[String],
) -> Result<Selected> {
    let alleles = alleles_of(indexes, calls, pools, contig_lengths, locus_ids);
    let joined = joined_allele(indexes, calls, pools, locus_ids);
    let mut with_join = alleles.clone();
    let annotated: FxSet<_> = indexes
        .iter()
        .filter_map(|&index| (!calls[index].refound).then_some(calls[index].dna))
        .collect();
    let longest_annotated = annotated
        .iter()
        .map(|&dna| &pools.dna[dna as usize])
        .filter(|dna| is_intact(dna))
        .map(String::len)
        .max()
        .unwrap_or(0);
    let mut excluded_refound_short = 0;
    if longest_annotated > 0 {
        let annotated_sequences: FxSet<_> = annotated
            .iter()
            .map(|&dna| pools.dna[dna as usize].as_str())
            .collect();
        with_join.retain(|allele| {
            let keep = annotated_sequences.contains(allele.dna.as_str())
                || allele.dna.len() * 100 >= longest_annotated * 95;
            if !keep {
                excluded_refound_short += 1;
                log::debug!(
                    "reference_candidate\texcluded_refound_short\tcluster\t{cluster}\tbp\t{}",
                    allele.dna.len()
                );
            }
            keep
        });
    }
    if let Some(dna) = joined {
        let longer = alleles.iter().all(|allele| dna.len() > allele.dna.len());
        if longer {
            with_join.push(Allele {
                dna,
                count: 0,
                synthetic: true,
                finished: 0,
            });
        }
    }
    let choice = choose_reference(&with_join, others, cluster, contig_lengths.len())
        .with_context(|| format!("No reference candidate for final family {cluster}"))?;
    let exact_samples: HashSet<_> = indexes
        .iter()
        .filter_map(|&index| {
            (pools.dna[calls[index].dna as usize] == choice.dna).then_some(calls[index].sample)
        })
        .collect();
    let mut sources: Vec<_> = indexes
        .iter()
        .filter_map(|&index| {
            (pools.dna[calls[index].dna as usize] == choice.dna).then(|| member_ids[index].clone())
        })
        .collect();
    if choice.selection == "joined_fragments" {
        for (dna, parts) in joined_models(indexes, calls, pools, locus_ids) {
            if dna == choice.dna {
                sources.extend(parts.into_iter().map(|index| member_ids[index].clone()));
            }
        }
    }
    sources.sort();
    Ok(Selected {
        excluded_refound_short,
        distinct: alleles.len(),
        intact: alleles
            .iter()
            .filter(|allele| is_intact(&allele.dna))
            .count(),
        exact_samples: exact_samples.len(),
        choice,
        sources,
    })
}

fn read_members(
    path: &str,
    families: &FxMap<&str, (u32, u32)>,
    sample_count: usize,
) -> Result<Vec<Member>> {
    let mut members = Vec::with_capacity(families.len());
    let mut seen = vec![false; families.len()];
    let mut seen_count = 0usize;
    let mut lines = LineReader::open(path)?;
    let header = lines.next_line().context("Empty reference metadata")??;
    if !header.starts_with(
        "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,",
    ) {
        bail!("Reference selection requires the eight-column gene_data.csv schema");
    }
    let mut line_number = 0usize;
    while let Some(line) = lines.next_line() {
        let line = line?;
        let mut fields = [""; 8];
        let mut count = 0;
        for field in line.splitn(8, ',') {
            fields[count] = field;
            count += 1;
        }
        let current_line = line_number;
        line_number += 1;
        if count < 3 {
            continue;
        }
        let Some(&(cluster, slot)) = families.get(fields[2]) else {
            continue;
        };
        if count != 8 {
            bail!("Invalid reference metadata at {}:{}", path, current_line + 2);
        }
        let id = fields[2];
        if std::mem::replace(&mut seen[slot as usize], true) {
            bail!("Duplicate reference member ID {id}");
        }
        seen_count += 1;
        let sample: usize = id
            .split('_')
            .next()
            .context("Missing member genome ID")?
            .parse()
            .with_context(|| format!("Invalid genome ID in {id}"))?;
        if sample >= sample_count {
            bail!("Genome ID out of range in reference member {id}");
        }
        let dna = uppercase_dna(fields[5]);
        if dna.is_empty() {
            bail!("Empty reference member sequence {id}");
        }
        let refound = id.contains("refound");
        members.push(Member {
            id: id.into(),
            sample,
            cluster,
            scaffold: fields[1].into(),
            annotation: fields[3].into(),
            dna,
            refound,
            location: if refound {
                Some(refound_location(fields[7], id)?)
            } else {
                None
            },
        });
    }
    if seen_count != families.len() {
        if let Some(missing) = families
            .iter()
            .find_map(|(id, &(_, slot))| (!seen[slot as usize]).then_some(*id))
        {
            bail!("Missing metadata for reference member {missing}");
        }
    }
    Ok(members)
}

fn refound_location(description: &str, id: &str) -> Result<Location> {
    let (coordinates, strand) = description
        .strip_prefix("location:")
        .and_then(|text| text.split_once(";strand:"))
        .with_context(|| format!("Missing refound coordinates for {id}"))?;
    let (start, end) = coordinates
        .split_once(':')
        .or_else(|| coordinates.split_once('-'))
        .with_context(|| format!("Invalid refound coordinates for {id}"))?;
    let start: u32 = start
        .parse()
        .with_context(|| format!("Invalid refound start for {id}"))?;
    let end: u32 = end
        .parse()
        .with_context(|| format!("Invalid refound end for {id}"))?;
    let strand = match strand {
        "+" => 1,
        "-" => -1,
        _ => bail!("Invalid refound strand for {id}"),
    };
    // Rustaroo stores refound locations as zero-based, half-open intervals.
    Ok(Location {
        start: start
            .checked_add(1)
            .context("Refound coordinate overflow")?,
        end,
        strand,
    })
}

fn read_locations(path: &str, needed: &NeededAnnotations) -> Result<SampleLocations> {
    let mut contigs = FxMap::default();
    let mut lengths = Vec::<usize>::new();
    let mut annotations = FxMap::default();
    let mut fasta = false;
    let mut current = None;
    let mut lines = LineReader::open(path)?;
    let mut number = 0usize;
    while let Some(line) = lines.next_line() {
        let line = line?.trim();
        let current_number = number;
        number += 1;
        if line == "##FASTA" {
            if fasta {
                bail!("Multiple FASTA sections in reference input {path}");
            }
            fasta = true;
            continue;
        }
        if !fasta {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            // Match the existing Rustaroo parser's removal of commas, but never
            // use this evidence parser to change accepted CDS or their sequence.
            let line: std::borrow::Cow<str> = if line.contains(',') {
                line.replace(',', "").into()
            } else {
                line.into()
            };
            // The first nine tab-separated fields, as a full split would give.
            let mut fields = [""; 9];
            let mut count = 0;
            for field in line.split('\t').take(9) {
                fields[count] = field;
                count += 1;
            }
            if count < 9 || fields[2] != "CDS" {
                continue;
            }
            let Some(id) = fields[8]
                .split(';')
                .find_map(|attr| attr.trim().strip_prefix("ID="))
            else {
                continue;
            };
            if !needed
                .get(fields[0])
                .is_some_and(|annotations| annotations.contains(id))
            {
                continue;
            }
            let key = (fields[0].to_string(), id.to_string());
            let number = current_number;
            let location = Location {
                start: fields[3]
                    .parse()
                    .with_context(|| format!("Invalid GFF start in {path}:{}", number + 1))?,
                end: fields[4]
                    .parse()
                    .with_context(|| format!("Invalid GFF end in {path}:{}", number + 1))?,
                strand: if fields[6] == "-" { -1 } else { 1 },
            };
            if annotations.insert(key, location).is_some() {
                bail!("Duplicate reference annotation {id} in {path}");
            }
        } else if let Some(header) = line.strip_prefix('>') {
            let id = header
                .split_whitespace()
                .next()
                .context("Empty reference contig name")?;
            if contigs.contains_key(id) {
                bail!("Duplicate reference contig {id} in {path}");
            }
            let index = u32::try_from(lengths.len())?;
            contigs.insert(id.to_string(), index);
            lengths.push(0);
            current = Some(index as usize);
        } else if !line.is_empty() {
            let index =
                current.with_context(|| format!("Sequence without a FASTA header in {path}"))?;
            lengths[index] = lengths[index]
                .checked_add(line.len())
                .context("Contig length overflow")?;
        }
    }
    if !fasta {
        bail!("Missing FASTA section in reference input {path}");
    }
    Ok(SampleLocations {
        contigs,
        lengths,
        annotations,
    })
}
