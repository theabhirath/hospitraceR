# Package index

## Transmission cluster detection algorithms

A collection of algorithms for detecting transmission clusters in
epidemiological data.

- [`get_tn_clusters_snp_thresh()`](https://theabhirath.github.io/hospitraceR/reference/get_tn_clusters_snp_thresh.md)
  : Cluster isolates using a hard SNP distance cutoff
- [`get_tn_clusters_sv_index()`](https://theabhirath.github.io/hospitraceR/reference/get_tn_clusters_sv_index.md)
  : Identify transmission clusters based on the number of shared
  variants

## Statistics

Functions for statistical analysis of transmission clustering results.

- [`cluster_pairwise_distances()`](https://theabhirath.github.io/hospitraceR/reference/cluster_pairwise_distances.md)
  : Compute intra-cluster SNP-distance summaries
- [`cluster_inter_distances()`](https://theabhirath.github.io/hospitraceR/reference/cluster_inter_distances.md)
  : Compute per-cluster nearest-neighbor SNP distances
- [`compute_genetic_distances_by_cluster()`](https://theabhirath.github.io/hospitraceR/reference/compute_genetic_distances_by_cluster.md)
  : Compute genetic distances by cluster
- [`intra_cluster_duration_metrics()`](https://theabhirath.github.io/hospitraceR/reference/intra_cluster_duration_metrics.md)
  : Compute intra-cluster duration metrics
- [`compute_duration_metrics_by_cluster()`](https://theabhirath.github.io/hospitraceR/reference/compute_duration_metrics_by_cluster.md)
  : Compute duration metrics by cluster

## Cluster comparison

Functions for comparing different cluster assignments using contingency
tables, information-theoretic measures, and source concordance.

- [`cluster_contingency_table()`](https://theabhirath.github.io/hospitraceR/reference/cluster_contingency_table.md)
  : Create a contingency table from two cluster assignments
- [`adjusted_mutual_information()`](https://theabhirath.github.io/hospitraceR/reference/adjusted_mutual_information.md)
  : Compute the adjusted mutual information from cluster assignments
- [`ami_from_contingency()`](https://theabhirath.github.io/hospitraceR/reference/ami_from_contingency.md)
  : Compute the adjusted mutual information from a contingency table
- [`adjusted_rand_index()`](https://theabhirath.github.io/hospitraceR/reference/adjusted_rand_index.md)
  : Compute the adjusted Rand index from cluster assignments
- [`ari_from_contingency()`](https://theabhirath.github.io/hospitraceR/reference/ari_from_contingency.md)
  : Compute the adjusted Rand index from a contingency table
- [`fraction_convert_same_source()`](https://theabhirath.github.io/hospitraceR/reference/fraction_convert_same_source.md)
  : Compute the fraction of converts with the same source across cluster
  assignments
- [`fraction_convert_same_source_from_lookups()`](https://theabhirath.github.io/hospitraceR/reference/fraction_convert_same_source_from_lookups.md)
  : Compute the fraction of converts with the same source from isolate
  lookups

## Epidemiological overlap analysis

Functions for analyzing patient overlap within and between transmission
clusters, and categorizing clusters by epidemiological patterns.

- [`isolate_isolate_overlap()`](https://theabhirath.github.io/hospitraceR/reference/isolate_isolate_overlap.md)
  : Calculate isolate-isolate co-location overlap
- [`isolate_isolate_sequential_overlap()`](https://theabhirath.github.io/hospitraceR/reference/isolate_isolate_sequential_overlap.md)
  : Calculate isolate-isolate sequential overlap
- [`cluster_isolate_overlap()`](https://theabhirath.github.io/hospitraceR/reference/cluster_isolate_overlap.md)
  : Summarize isolate overlap within each cluster
- [`fraction_convert_events_with_overlap()`](https://theabhirath.github.io/hospitraceR/reference/fraction_convert_events_with_overlap.md)
  : Calculate the fraction of convert events with overlap
- [`categorize_cluster_overlap()`](https://theabhirath.github.io/hospitraceR/reference/categorize_cluster_overlap.md)
  : Categorize cluster overlap
- [`flatten_cluster_overlap_categorization()`](https://theabhirath.github.io/hospitraceR/reference/flatten_cluster_overlap_categorization.md)
  : Flatten cluster-overlap categorization to a data frame
- [`cluster_overlap_perm_test()`](https://theabhirath.github.io/hospitraceR/reference/cluster_overlap_perm_test.md)
  : Permutation test for cluster overlap fractions
- [`cluster_patient_categorization()`](https://theabhirath.github.io/hospitraceR/reference/cluster_patient_categorization.md)
  : Categorize patients within clusters by epidemiological status
- [`flatten_cluster_patient_categorization()`](https://theabhirath.github.io/hospitraceR/reference/flatten_cluster_patient_categorization.md)
  : Flatten cluster-patient categorization to a data frame

## Utilities

A collection of utility functions for this package.

- [`get_snp_dist_matrix()`](https://theabhirath.github.io/hospitraceR/reference/get_snp_dist_matrix.md)
  : Compute a pairwise SNP distance matrix from a DNA alignment
- [`get_phylo_tree()`](https://theabhirath.github.io/hospitraceR/reference/get_phylo_tree.md)
  : Build a phylogenetic tree by neighbor-joining or maximum parsimony
- [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md)
  : Create a lookup table for isolates along with their cluster
  assignments
- [`get_non_single_patient_clusters()`](https://theabhirath.github.io/hospitraceR/reference/get_non_single_patient_clusters.md)
  : Get clusters containing more than one patient
- [`get_singleton_acq_clusters()`](https://theabhirath.github.io/hospitraceR/reference/get_singleton_acq_clusters.md)
  : Get single-patient clusters that represent an acquisition
- [`remove_singleton_clusters()`](https://theabhirath.github.io/hospitraceR/reference/remove_singleton_clusters.md)
  : Remove singleton clusters from a vector of cluster assignments
