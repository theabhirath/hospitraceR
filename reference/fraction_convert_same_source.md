# Compute the fraction of converts with the same source across cluster assignments

Compares two cluster assignments by examining patients categorized as
"convert". For each patient who is a convert in both assignments, checks
whether the source patient (the cluster's index, weak index, or
multiply-colonized index) is the same across both assignments.

## Usage

``` r
fraction_convert_same_source(
  clusters1,
  clusters2,
  dna_aln,
  seq2pt,
  adm_seqs,
  dates,
  surv_df,
  converts_without_assigned_source = FALSE,
  debug = FALSE
)
```

## Arguments

- clusters1:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

- clusters2:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

- dna_aln:

  A DNA alignment object of class `DNAbin`.

- seq2pt:

  A named vector mapping sequence IDs to patient IDs.

- adm_seqs:

  A vector of sequence IDs corresponding to admission positive patients.

- dates:

  A named vector mapping sequence IDs to dates.

- surv_df:

  A data frame with surveillance data containing columns: patient_id,
  genome_id, surv_date, result.

- converts_without_assigned_source:

  A logical value indicating whether to include singleton converts –
  patients counted as converts that have no assigned source in the
  lookup – in the comparison, recording them with no source.

## Value

A numeric value between 0 and 1 representing the fraction of convert
patients whose source (index, weak index, or multiply-colonized index)
is the same in both cluster assignments. Returns NA if there are no
common converts.

## Details

The function:

1.  Creates isolate lookups for both cluster assignments

2.  Categorizes patients in both using
    [`cluster_patient_categorization()`](https://theabhirath.github.io/hospitraceR/reference/cluster_patient_categorization.md)

3.  Identifies patients categorized as "convert" in both assignments

4.  For each such convert, finds the source patient (index, weak index,
    or multiply-colonized index) in their cluster for each assignment

5.  Returns the fraction where the source patient matches

## See also

[`fraction_convert_same_source_from_lookups()`](https://theabhirath.github.io/hospitraceR/reference/fraction_convert_same_source_from_lookups.md),
[`cluster_patient_categorization()`](https://theabhirath.github.io/hospitraceR/reference/cluster_patient_categorization.md),
[`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md)
