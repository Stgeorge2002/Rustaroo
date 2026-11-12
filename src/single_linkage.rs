//! Single linkage clustering on a distance matrix

use petgraph::graph::NodeIndex;

pub fn single_linkage_clustering(
    distance_matrix: &[Vec<f64>],
    nodes: &[NodeIndex],
    threshold: f64,
) -> Vec<Vec<NodeIndex>> {
    if distance_matrix.is_empty() {
        return vec![];
    }

    single_linkage_clustering_by(nodes, |first, second| {
        distance_matrix
            .get(first)
            .and_then(|row| row.get(second))
            .is_some_and(|&distance| distance < threshold)
    })
}

/// Cluster nodes from a lazy pairwise predicate without materializing a matrix.
pub fn single_linkage_clustering_by<F>(
    nodes: &[NodeIndex],
    mut pair_is_close: F,
) -> Vec<Vec<NodeIndex>>
where
    F: FnMut(usize, usize) -> bool,
{
    if nodes.is_empty() {
        return vec![];
    }

    let mut parent: Vec<usize> = (0..nodes.len()).collect();
    let mut rank = vec![0usize; nodes.len()];

    for i in 0..nodes.len() {
        for j in (i + 1)..nodes.len() {
            if find(&mut parent, i) != find(&mut parent, j) && pair_is_close(i, j) {
                union(&mut parent, &mut rank, i, j);
            }
        }
    }

    let mut root_to_cluster = vec![usize::MAX; nodes.len()];
    let mut clusters: Vec<Vec<NodeIndex>> = Vec::new();
    for (index, &node) in nodes.iter().enumerate() {
        let root = find(&mut parent, index);
        let cluster_index = if root_to_cluster[root] == usize::MAX {
            let cluster_index = clusters.len();
            root_to_cluster[root] = cluster_index;
            clusters.push(Vec::new());
            cluster_index
        } else {
            root_to_cluster[root]
        };
        clusters[cluster_index].push(node);
    }

    clusters
}

fn find(parent: &mut [usize], mut item: usize) -> usize {
    let mut root = item;
    while parent[root] != root {
        root = parent[root];
    }
    while parent[item] != item {
        let next = parent[item];
        parent[item] = root;
        item = next;
    }
    root
}

fn union(parent: &mut [usize], rank: &mut [usize], first: usize, second: usize) {
    let mut first_root = find(parent, first);
    let mut second_root = find(parent, second);
    if first_root == second_root {
        return;
    }
    if rank[first_root] < rank[second_root] {
        std::mem::swap(&mut first_root, &mut second_root);
    }
    parent[second_root] = first_root;
    if rank[first_root] == rank[second_root] {
        rank[first_root] += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reference_clustering(
        distance_matrix: &[Vec<f64>],
        nodes: &[NodeIndex],
        threshold: f64,
    ) -> Vec<Vec<NodeIndex>> {
        if nodes.is_empty() || distance_matrix.is_empty() {
            return vec![];
        }
        let mut clusters: Vec<Vec<usize>> = (0..nodes.len()).map(|i| vec![i]).collect();
        loop {
            let mut best_pair = None;
            let mut best_distance = threshold;
            for i in 0..clusters.len() {
                for j in (i + 1)..clusters.len() {
                    let mut minimum = f64::MAX;
                    for &first in &clusters[i] {
                        for &second in &clusters[j] {
                            if let Some(&distance) =
                                distance_matrix.get(first).and_then(|row| row.get(second))
                            {
                                minimum = minimum.min(distance);
                            }
                        }
                    }
                    if minimum < best_distance {
                        best_distance = minimum;
                        best_pair = Some((i, j));
                    }
                }
            }
            match best_pair {
                Some((i, j)) => {
                    let merged = clusters.remove(j);
                    clusters[i].extend(merged);
                }
                None => break,
            }
        }
        clusters
            .into_iter()
            .map(|cluster| cluster.into_iter().map(|index| nodes[index]).collect())
            .collect()
    }

    fn canonicalize(mut clusters: Vec<Vec<NodeIndex>>) -> Vec<Vec<NodeIndex>> {
        for cluster in &mut clusters {
            cluster.sort();
        }
        clusters.sort_by_key(|cluster| cluster.first().copied());
        clusters
    }

    #[test]
    fn union_find_matches_iterative_single_linkage() {
        let nodes: Vec<_> = (0..5).map(NodeIndex::new).collect();
        let matrices = [
            vec![
                vec![0.0, 0.1, 0.9, 0.9, 0.9],
                vec![0.1, 0.0, 0.2, 0.9, 0.9],
                vec![0.9, 0.2, 0.0, 0.8, 0.9],
                vec![0.9, 0.9, 0.8, 0.0, 0.3],
                vec![0.9, 0.9, 0.9, 0.3, 0.0],
            ],
            vec![
                vec![0.0, 0.5, 0.1, 0.9, 0.9],
                vec![0.5, 0.0, 0.4, 0.1, 0.9],
                vec![0.1, 0.4, 0.0, 0.9, 0.9],
                vec![0.9, 0.1, 0.9, 0.0, 0.2],
                vec![0.9, 0.9, 0.9, 0.2, 0.0],
            ],
        ];

        for matrix in matrices {
            for threshold in [0.1, 0.3, 0.5, 0.8, 1.0] {
                assert_eq!(
                    canonicalize(single_linkage_clustering(&matrix, &nodes, threshold)),
                    canonicalize(reference_clustering(&matrix, &nodes, threshold))
                );
            }
        }
    }
}
