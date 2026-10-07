# Categorize patients within clusters by epidemiological status

Categorizes each patient in a cluster based on their role in
transmission:

- index: admission positive, earliest isolate, isolate in cluster

- multiply-colonized-index: admission positive, earliest isolate,
  isolate not in cluster

- weak-index: not admission positive, but first surveillance is positive
  and in cluster

- convert: had surveillance before first positive

- adm-pos: first positive is in cluster, is first surveillance, and is
  admission positive

- adm-pos-convert: first positive is not in cluster but is first
  surveillance

- secondary-convert: first positive is not in cluster and had prior
  surveillance

- ambiguous-adm-pos: first positive is culture-only (strain unknown) and
  is first surveillance

- ambiguous-convert: first positive is culture-only (strain unknown) and
  had prior surveillance, or is in cluster and is first surveillance but
  not admission positive

## Usage

``` r
cluster_patient_categorization(isolate_lookup, surv_df)
```

## Arguments

- isolate_lookup:

  A lookup table from
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md).

- surv_df:

  A data frame with surveillance data (patient_id, genome_id, surv_date,
  result).

## Value

A named list where each element is a named character vector mapping
patient_id to their category within that cluster.
