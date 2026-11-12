//! Tests for network module

use rustaroo::network::GeneNode;

#[test]
fn test_gene_node_creation() {
    let node = GeneNode::default();
    assert_eq!(node.size, 0);
    assert!(!node.paralog);
    assert!(!node.has_end);
    assert!(node.members.is_empty());
    assert!(node.seq_ids.is_empty());
}

#[test]
fn test_gene_node_with_data() {
    let mut node = GeneNode::default();
    node.size = 5;
    node.annotation = "test_gene".to_string();
    node.paralog = true;

    assert_eq!(node.size, 5);
    assert_eq!(node.annotation, "test_gene");
    assert!(node.paralog);
}
