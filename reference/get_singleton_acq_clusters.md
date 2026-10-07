# Get single-patient clusters that represent an acquisition

The complement of
[`get_non_single_patient_clusters()`](https://theabhirath.github.io/hospitraceR/reference/get_non_single_patient_clusters.md),
restricted to acquisitions: clusters made up of isolates from exactly
one patient, where that patient has no admission-positive isolate.
Singleton clusters of admission-positive patients are excluded, since
they were not acquired in the facility and so nothing about them needs
an overlap explanation.

## Usage

``` r
get_singleton_acq_clusters(isolate_lookup)
```

## Arguments

- isolate_lookup:

  A lookup table for isolates and their cluster assignments with other
  relevant epidemiological information. See
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md).

## Value

A numeric vector of cluster IDs containing isolates from exactly one,
non admission-positive, patient.
