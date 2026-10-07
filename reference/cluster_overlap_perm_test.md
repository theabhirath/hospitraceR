# Permutation test for cluster overlap fractions

Performs a permutation test to assess whether the observed fraction of
converts with overlap is significantly different from what would be
expected by chance. Tests overlap at facility, floor, and room levels.

## Usage

``` r
cluster_overlap_perm_test(
  clusters,
  dna_aln,
  seq2pt,
  adm_seqs,
  adm_pos_pt_seqs,
  dates,
  surv_df,
  facility_trace,
  floor_trace,
  room_trace,
  nperm = 1000,
  num_cores = detectCores() - 1,
  include_singleton_acq = FALSE
)
```

## Arguments

- clusters:

  A named numeric vector of cluster assignments.

- dna_aln:

  A DNA alignment object (used for isolate IDs).

- seq2pt:

  A named vector mapping sequence IDs to patient IDs.

- adm_seqs:

  A vector of sequence IDs which correspond to admission positive
  sequences.

- adm_pos_pt_seqs:

  A vector of sequence IDs for patients who are admission positive.

- dates:

  A named vector mapping sequence IDs to dates.

- surv_df:

  A data frame with surveillance data.

- facility_trace:

  A matrix with patient IDs as row names and dates as column names.

- floor_trace:

  A matrix with floor-level trace data (same structure as
  `facility_trace`).

- room_trace:

  A matrix with room-level trace data (same structure as
  `facility_trace`).

- nperm:

  Number of permutations to perform.

- num_cores:

  Number of cores for parallel processing.

- include_singleton_acq:

  Whether to evaluate all acquisitions, including those not classified
  to be a part of multi-patient clusters.

## Value

A list containing:

- `observed`: A matrix with one row per trace type (facility, floor,
  room, seq_facility, seq_floor, seq_room) and columns `n_overlap`
  (converts with overlap) and `n_converts` (converts), each summed over
  the tested clusters

- `permuted`: A numeric array of dimensions (trace type, count, nperm)
  holding the same two sums for every permutation

- `valid_clusters`: A numeric vector of the cluster IDs tested — those
  with more than one patient, plus the single-patient acquisition
  clusters when `include_singleton_acq` is TRUE

The pooled, convert-weighted overlap fraction is
`n_overlap / n_converts`; summing the counts across sequence types
before dividing pools them further.

A permutation whose cluster assignment strands (capacity too
concentrated to place every patient) stops the whole test with an error.
