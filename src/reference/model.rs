//! Reference evidence only. Cluster IDs are assigned by the existing Rustaroo graph.

pub(super) struct Call {
    pub sample: u32,
    pub contig: u32,
    pub start: u32,
    pub end: u32,
    pub strand: i8,
    pub dna: u32,
    pub protein: u32,
    pub cluster: u32,
    pub locus: Option<u32>,
    pub refound: bool,
}

#[derive(Default)]
pub(super) struct Pools {
    pub dna: Vec<String>,
    pub protein: Vec<Vec<u8>>,
}
