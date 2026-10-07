# Identify transmission clusters based on the number of shared variants

Clustering is performed to identify the maximal clusters containing a
single intake-positive patient that occurs before all cluster converts.
The clustering metric is the number of shared variants, and clusters can
have multiple intake-positive patients if they share an identical number
of variants with other cluster members or intake-positive patients occur
after converts. This clustering also requires that clusters be defined
by at least one shared variant that other isolates don't have.

## Usage

``` r
get_tn_clusters_sv_index(
  dna_aln,
  snp_dist,
  adm_seqs,
  adm_pos_pt_seqs,
  seq2pt,
  dates,
  tree,
  snp_thresh = Inf
)
```

## Arguments

- dna_aln:

  A DNA alignment object of class `DNAbin`. Its first sequence is taken
  as the outgroup.

- snp_dist:

  A matrix of SNP distances between isolates. See
  [`get_snp_dist_matrix()`](https://theabhirath.github.io/hospitraceR/reference/get_snp_dist_matrix.md).

- adm_seqs:

  A vector of sequence IDs for sequences from patients positive at
  intake.

- adm_pos_pt_seqs:

  A vector of all sequence IDs from admission-positive patients, whether
  collected at intake or later; a superset of `adm_seqs`.

- seq2pt:

  A named vector mapping sequence IDs to patient IDs.

- dates:

  A named vector of isolate dates, named by sequence ID.

- tree:

  A phylogenetic tree of class `phylo` over the same isolates, e.g. from
  [`get_phylo_tree()`](https://theabhirath.github.io/hospitraceR/reference/get_phylo_tree.md).
  Its first tip is taken as the outgroup, matching `dna_aln`.

- snp_thresh:

  SNP distance above which a threshold-free cluster is broken up. Within
  each cluster, isolates connected through a chain of pairwise distances
  at or below this value (single-linkage components) are kept together;
  disconnected groups are split into separate clusters, all retained.
  `Inf` (the default) never breaks a cluster.

## Value

A numeric vector giving the cluster each isolate belongs to.

## References

Hawken, S. E., Yelin, R. D., Lolans, K., Pirani, A., Weinstein, R. A.,
Lin, M. Y., Hayden, M. K., & Snitkin, E. S. (2022). Threshold-free
genomic cluster detection to track transmission pathways in health-care
settings: A genomic epidemiology analysis. The Lancet Microbe, 3(9),
e652–e662.
[doi:10.1016/S2666-5247(22)00115-X](https://doi.org/10.1016/S2666-5247%2822%2900115-X)
