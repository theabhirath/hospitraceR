# Compute the fraction of converts with the same source from isolate lookups

Compares two isolate lookups by examining patients categorized as
"convert". This works directly with isolate lookups and is useful when
the two cluster assignments do not share the same epidemiological
information.

## Usage

``` r
fraction_convert_same_source_from_lookups(
  isolate_lookup1,
  isolate_lookup2,
  surv_df_1,
  surv_df_2,
  converts_without_assigned_source = TRUE,
  debug = FALSE
)
```

## Arguments

- isolate_lookup1:

  A data frame from
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md)
  for the first cluster assignment.

- isolate_lookup2:

  A data frame from
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md)
  for the second cluster assignment.

- surv_df_1:

  A data frame with surveillance data containing columns: patient_id,
  genome_id, surv_date, result for the first cluster assignment.

- surv_df_2:

  A data frame with surveillance data containing columns: patient_id,
  genome_id, surv_date, result for the second cluster assignment.

- converts_without_assigned_source:

  A logical value indicating whether to include singleton converts –
  patients counted as converts that have no assigned source in the
  lookup – in the comparison, recording them with no source.

## Value

A numeric value between 0 and 1 representing the fraction of convert
patients whose source (index, weak index, or multiply-colonized index)
is the same in both cluster assignments. Returns NA if there are no
common converts.

## See also

[`fraction_convert_same_source()`](https://theabhirath.github.io/hospitraceR/reference/fraction_convert_same_source.md),
[`cluster_patient_categorization()`](https://theabhirath.github.io/hospitraceR/reference/cluster_patient_categorization.md)
