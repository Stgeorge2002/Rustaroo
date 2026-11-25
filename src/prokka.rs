//! Complete Prokka GFF parsing and sequence extraction
//!

use anyhow::{Context, Result};
use bio::io::fasta;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Gene record extracted from GFF
#[derive(Debug, Clone)]
pub struct GeneRecord {
    pub id: String,
    pub annotation_id: String,
    pub name: String,
    pub description: String,
    pub scaffold: String,
    pub start: usize,
    pub dna_sequence: String,
    pub protein_sequence: String,
}

/// Translate DNA to protein using bacterial codon table
pub fn translate(dna: &str, table: usize) -> String {
    // ASCII input (all real DNA) takes a byte-level path with exactly the same
    // result; anything else uses the original Unicode-aware implementation.
    if !dna.is_ascii() {
        return translate_str(dna, table);
    }
    let dna = dna.as_bytes();
    let mut protein = Vec::with_capacity(dna.len() / 3);
    for codon in dna.chunks_exact(3) {
        protein.push(codon_amino([
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ]));
    }
    if table == 11
        && dna.len() >= 3
        && matches!(
            [
                dna[0].to_ascii_uppercase(),
                dna[1].to_ascii_uppercase(),
                dna[2].to_ascii_uppercase(),
            ],
            [b'T', b'T', b'G']
                | [b'C', b'T', b'G']
                | [b'A', b'T', b'T']
                | [b'A', b'T', b'C']
                | [b'A', b'T', b'A']
                | [b'A', b'T', b'G']
                | [b'G', b'T', b'G']
        )
        && !protein.is_empty()
    {
        protein[0] = b'M';
    }
    String::from_utf8(protein).expect("translation is ASCII")
}

#[inline]
fn codon_amino(codon: [u8; 3]) -> u8 {
    match &codon {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X', // Unknown codon
    }
}

/// Reference translation for any (including non-ASCII) input.
fn translate_str(dna: &str, table: usize) -> String {
    let dna_upper = dna.to_uppercase();
    let mut protein = String::new();

    for codon_start in (0..dna_upper.len()).step_by(3) {
        if codon_start + 3 > dna_upper.len() {
            break;
        }

        let codon = &dna_upper[codon_start..codon_start + 3];
        let aa = match codon {
            // Use standard bacterial translation
            "TTT" | "TTC" => "F",
            "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => "L",
            "ATT" | "ATC" | "ATA" => "I",
            "ATG" => "M",
            "GTT" | "GTC" | "GTA" | "GTG" => "V",
            "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => "S",
            "CCT" | "CCC" | "CCA" | "CCG" => "P",
            "ACT" | "ACC" | "ACA" | "ACG" => "T",
            "GCT" | "GCC" | "GCA" | "GCG" => "A",
            "TAT" | "TAC" => "Y",
            "TAA" | "TAG" | "TGA" => "*",
            "CAT" | "CAC" => "H",
            "CAA" | "CAG" => "Q",
            "AAT" | "AAC" => "N",
            "AAA" | "AAG" => "K",
            "GAT" | "GAC" => "D",
            "GAA" | "GAG" => "E",
            "TGT" | "TGC" => "C",
            "TGG" => "W",
            "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => "R",
            "GGT" | "GGC" | "GGA" | "GGG" => "G",
            _ => "X", // Unknown codon
        };
        protein.push_str(aa);
    }

    if table == 11
        && dna_upper.len() >= 3
        && matches!(
            &dna_upper[..3],
            "TTG" | "CTG" | "ATT" | "ATC" | "ATA" | "ATG" | "GTG"
        )
        && !protein.is_empty()
    {
        protein.replace_range(..1, "M");
    }
    protein
}

/// Reverse complement a DNA sequence
pub fn reverse_complement(dna: &str) -> String {
    if !dna.is_ascii() {
        return reverse_complement_chars(dna);
    }
    // Byte-level fast path; identical to the char-level version for ASCII.
    let bytes: Vec<u8> = dna
        .bytes()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            _ => c,
        })
        .collect();
    String::from_utf8(bytes).expect("reverse complement of ASCII is ASCII")
}

fn reverse_complement_chars(dna: &str) -> String {
    dna.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            _ => c,
        })
        .collect()
}

/// Parse GFF3 and extract gene sequences
pub fn get_gene_sequences(
    gff_file_name: &str,
    file_number: usize,
    filter_seqs: bool,
    table: usize,
) -> Result<(HashMap<String, GeneRecord>, Vec<GeneRecord>)> {
    let protein_list = get_gene_list(gff_file_name, file_number, filter_seqs, table)?;
    let sequence_dictionary = protein_list
        .iter()
        .map(|record| (record.id.clone(), record.clone()))
        .collect();
    Ok((sequence_dictionary, protein_list))
}

/// Parse GFF3 and extract gene records in clustering-id order, without the
/// per-id dictionary built by [`get_gene_sequences`].
pub fn get_gene_list(
    gff_file_name: &str,
    file_number: usize,
    filter_seqs: bool,
    table: usize,
) -> Result<Vec<GeneRecord>> {
    let content = std::fs::read_to_string(gff_file_name)
        .with_context(|| format!("Failed to read GFF file: {}", gff_file_name))?;

    let content = content.replace(",", "");
    let parts: Vec<&str> = content.split("##FASTA").collect();

    if parts.len() != 2 {
        anyhow::bail!("GFF file must have exactly one ##FASTA section");
    }

    let gff_section = parts[0];
    let fasta_section = parts[1];

    // Parse FASTA sequences
    let mut sequences = HashMap::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in fasta_section.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if !current_id.is_empty() {
                sequences.insert(
                    std::mem::take(&mut current_id),
                    std::mem::take(&mut current_seq),
                );
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
        sequences.insert(current_id, current_seq);
    }

    // Parse GFF entries
    let mut scaffold_genes: BTreeMap<String, Vec<(usize, GeneRecord)>> = BTreeMap::new();
    let mut scaffold_order = Vec::new();

    for line in gff_section.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let feature_type = fields[2];
        if feature_type != "CDS" {
            continue;
        }

        let seqid = fields[0];
        let start: usize = fields[3].parse().context("Invalid GFF start coordinate")?;
        let end: usize = fields[4].parse().context("Invalid GFF end coordinate")?;
        let strand = fields[6];
        let frame = fields[7];
        let attributes = fields[8];

        // Check frame
        if frame != "0" {
            if filter_seqs {
                continue;
            } else {
                anyhow::bail!("Frame shifts not supported");
            }
        }

        // Parse attributes
        let mut gene_id = String::new();
        let mut gene_name = String::new();
        let mut product = String::new();

        for attr in attributes.split(';') {
            let attr = attr.trim();
            if attr.starts_with("ID=") {
                gene_id = attr[3..].to_string();
            } else if attr.starts_with("gene=") {
                gene_name = attr[5..].to_string();
            } else if attr.starts_with("product=") {
                product = attr[8..].to_string();
            } else if attr.starts_with("name=") && gene_name.is_empty() {
                gene_name = attr[5..].to_string();
            }
        }

        // Get sequence
        let scaffold_seq = sequences.get(seqid);
        if scaffold_seq.is_none() {
            if filter_seqs {
                continue;
            } else {
                anyhow::bail!("Sequence ID not found: {}", seqid);
            }
        }

        let scaffold_seq = scaffold_seq.unwrap();
        if start == 0 || end < start || end > scaffold_seq.len() {
            anyhow::bail!(
                "GFF gene {} has out-of-range coordinates {}..{}",
                gene_id,
                start,
                end
            );
        }
        let mut gene_seq = scaffold_seq[start - 1..end].to_string();

        if strand == "-" {
            gene_seq = reverse_complement(&gene_seq);
        }

        // Validate sequence
        if gene_seq.len() % 3 != 0 || gene_seq.len() < 34 {
            if filter_seqs {
                continue;
            } else {
                anyhow::bail!("Invalid gene length: {}", gene_id);
            }
        }

        let protein_seq = translate(&gene_seq, table);

        // Check for premature stop codons
        if protein_seq[..protein_seq.len().saturating_sub(1)].contains('*') {
            if filter_seqs {
                continue;
            } else {
                anyhow::bail!("Premature stop codon in: {}", gene_id);
            }
        }

        let record = GeneRecord {
            id: gene_id.clone(),
            annotation_id: gene_id,
            name: gene_name,
            description: product,
            scaffold: seqid.to_string(),
            start,
            dna_sequence: gene_seq,
            protein_sequence: protein_seq.trim_end_matches('*').to_string(),
        };

        if !scaffold_genes.contains_key(seqid) {
            scaffold_order.push(seqid.to_string());
        }
        scaffold_genes
            .entry(seqid.to_string())
            .or_insert_with(Vec::new)
            .push((start, record));
    }

    // Sort genes by position and assign clustering IDs. Ids are unique
    // (file, scaffold, index), so the list has one record per id.
    let mut protein_list = Vec::new();

    let mut scaff_count = 0;
    for scaffold in scaffold_order {
        let mut genes = scaffold_genes.remove(&scaffold).unwrap();
        genes.sort_by_key(|(start, _)| *start);

        for (gene_index, (_start, mut record)) in genes.into_iter().enumerate() {
            record.id = format!("{}_{}_{}", file_number, scaff_count, gene_index);
            protein_list.push(record);
        }
        scaff_count += 1;
    }

    if protein_list.is_empty() {
        anyhow::bail!("No valid sequences found in GFF");
    }

    Ok(protein_list)
}

/// Write sequences to FASTA files and CSV
pub fn output_files(
    _dna_dict: &HashMap<String, GeneRecord>,
    protein_list: &[GeneRecord],
    prot_file: &str,
    dna_file: &str,
    csv_file: &str,
    gff_filename: &str,
) -> Result<()> {
    // Write protein sequences
    {
        let file = File::create(prot_file)?;
        let mut writer = fasta::Writer::new(file);
        for record in protein_list {
            let fasta_record =
                fasta::Record::with_attrs(&record.id, None, record.protein_sequence.as_bytes());
            writer.write_record(&fasta_record)?;
        }
    }

    // Write DNA sequences
    {
        let file = File::create(dna_file)?;
        let mut writer = fasta::Writer::new(file);
        for record in protein_list {
            let fasta_record =
                fasta::Record::with_attrs(&record.id, None, record.dna_sequence.as_bytes());
            writer.write_record(&fasta_record)?;
        }
    }

    // Write CSV metadata
    {
        let mut file = std::io::BufWriter::new(File::create(csv_file)?);
        writeln!(
            file,
            "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description"
        )?;

        let gff_name = Path::new(gff_filename)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("");

        for record in protein_list {
            writeln!(
                file,
                "{},{},{},{},{},{},{},{}",
                gff_name,
                record.scaffold,
                record.id,
                record.annotation_id,
                record.protein_sequence,
                record.dna_sequence,
                record.name,
                record.description
            )?;
        }
        file.flush()?;
    }

    Ok(())
}

/// Process all Prokka input files
pub fn process_prokka_input(
    input_files: &[String],
    filter_invalid: bool,
    codon_table: usize,
) -> Result<Vec<HashMap<String, GeneRecord>>> {
    let mut all_genes = Vec::new();

    for (i, input_file) in input_files.iter().enumerate() {
        let (genes, _proteins) = get_gene_sequences(input_file, i, filter_invalid, codon_table)?;
        all_genes.push(genes);
    }

    Ok(all_genes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate() {
        assert_eq!(translate("ATGGCCTAG", 11), "MA*");
        assert_eq!(translate("ATGAAATAG", 11), "MK*");
    }

    #[test]
    fn byte_paths_match_reference_implementations() {
        let alphabet = b"ACGTacgtNnRYXE*-";
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        for length in 0..200 {
            for _ in 0..20 {
                let dna: String = (0..length)
                    .map(|_| {
                        state ^= state << 13;
                        state ^= state >> 7;
                        state ^= state << 17;
                        alphabet[state as usize % alphabet.len()] as char
                    })
                    .collect();
                for table in [1, 4, 11] {
                    assert_eq!(translate(&dna, table), translate_str(&dna, table), "{dna}");
                }
                assert_eq!(reverse_complement(&dna), reverse_complement_chars(&dna));
            }
        }
        for dna in ["ATGß", "ÄTGAAA", "ttgaaataa", "ctg", "GTGAAATAG"] {
            assert_eq!(translate(dna, 11), translate_str(dna, 11));
            assert_eq!(reverse_complement(dna), reverse_complement_chars(dna));
        }
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("ATCG"), "CGAT");
    }
}
