# Categorize cluster overlap

Labels each cluster by its likely transmission route, from the earliest
collected isolate among all patients in the cluster (the "index
isolate") and whether the cluster's other isolates have an overlap
explanation. The categories are:

- "patient-to-patient": if the index isolate is admission-positive and
  there is overlap explanation for all other converts in the cluster.

- "weak-index-patient-to-patient": if the index isolate is not
  admission-positive but is the first surveillance for the patient after
  admission and there is overlap explanation for all other isolates in
  the cluster.

- "missing-intermediate": if the index isolate is admission-positive but
  at least one other convert in the cluster has no overlap explanation.

- "weak-index-missing-intermediate": if the index isolate is a weak
  index (not admission-positive but the first surveillance for the
  patient after admission) and at least one other convert in the cluster
  has no overlap explanation.

- "missing-source": if the index isolate is in the cluster but is
  neither admission-positive nor a weak index (the index isolate is not
  the patient's first surveillance culture), so the cluster's source was
  not sampled. The other isolates' overlap explanations are not
  consulted.

- "multiply-colonized-index": if the index isolate is not in the cluster
  but is admission-positive, this is a "multiply-colonized index" if
  there is overlap explanation for all other isolates in the cluster.

- "multiply-colonized-index-missing-intermediate": if the index isolate
  is not in the cluster but is admission-positive, and at least one
  other convert in the cluster has no overlap explanation.

- "all-admission-positive": if all isolates in the cluster are admission
  positive, this is a special category that may indicate a common source
  or multiple colonization events, but is not explained

- "inexplicable": catch-all category for cases that do not fit into the
  other categories. Currently, index and overlap explanations are not
  enough to explain these clusters.

## Usage

``` r
categorize_cluster_overlap(isolate_lookup, cluster_overlap_df, surv_df)
```

## Arguments

- isolate_lookup:

  A lookup table for isolates and their cluster assignments with other
  relevant epidemiological information. See
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md).

- cluster_overlap_df:

  A data frame of per-isolate cluster overlap, from
  [`cluster_isolate_overlap()`](https://theabhirath.github.io/hospitraceR/reference/cluster_isolate_overlap.md).

- surv_df:

  A data frame of surveillance cultures, with columns `patient_id`,
  `genome_id`, `surv_date` and `result`.

## Value

A named character vector giving the overlap category of each cluster.
