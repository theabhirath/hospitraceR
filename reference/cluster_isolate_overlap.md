# Summarize isolate overlap within each cluster

For each isolate in a cluster (treated as a recipient), reports whether
any other isolate in the same cluster (as a donor) overlaps with it,
using the isolate-pair overlaps from
[`isolate_isolate_overlap()`](https://theabhirath.github.io/hospitraceR/reference/isolate_isolate_overlap.md)
or
[`isolate_isolate_sequential_overlap()`](https://theabhirath.github.io/hospitraceR/reference/isolate_isolate_sequential_overlap.md).
Admission-positive isolates are reported as `NA`, since they were not
acquired in the facility.

Only multi-patient clusters are reported by default; set
`include_singleton_acq` to also account for acquisitions sitting in
single-patient clusters.

## Usage

``` r
cluster_isolate_overlap(
  isolate_lookup,
  iso_overlap_df,
  include_singleton_acq = FALSE
)
```

## Arguments

- isolate_lookup:

  A lookup table for isolates and their cluster assignments with other
  relevant epidemiological information. See
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md).

- iso_overlap_df:

  A data frame of isolate-pair overlaps, from
  [`isolate_isolate_overlap()`](https://theabhirath.github.io/hospitraceR/reference/isolate_isolate_overlap.md).

- include_singleton_acq:

  Whether to evaluate all acquisitions, including those not classified
  to be a part of multi-patient clusters.

## Value

A data frame with columns `cluster`, `isolate_id` and `overlap`.
