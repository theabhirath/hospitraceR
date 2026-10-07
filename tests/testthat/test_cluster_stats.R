# Tests for the per-cluster genetic-distance and duration metrics with synthetic inputs.

# A symmetric SNP distance matrix over named isolates, zero on the diagonal.
make_snp_dist <- function(labels, pairs) {
    m <- matrix(0, length(labels), length(labels), dimnames = list(labels, labels))
    for (nm in names(pairs)) {
        ij <- strsplit(nm, "-")[[1]]
        m[ij[1], ij[2]] <- pairs[[nm]]
        m[ij[2], ij[1]] <- pairs[[nm]]
    }
    m
}

# --- cluster_pairwise_distances ------------------------------------------------

test_that("cluster_pairwise_distances summarizes the upper triangle of a cluster", {
    snp_dist <- make_snp_dist(
        c("a", "b", "c"),
        list("a-b" = 2, "a-c" = 4, "b-c" = 6)
    )

    out <- cluster_pairwise_distances(c("a", "b", "c"), snp_dist)

    expect_equal(unname(out["mean_genetic_distance"]), 4) # (2 + 4 + 6) / 3
    expect_equal(unname(out["median_genetic_distance"]), 4)
    expect_equal(unname(out["max_genetic_distance"]), 6)
})

# --- cluster_inter_distances ---------------------------------------------------

test_that("cluster_inter_distances separates nearest other-cluster vs nearest other-isolate", {
    # a, b -> cluster 1; c, d -> cluster 2; x is an unclustered isolate present in the matrix.
    labels <- c("a", "b", "c", "d", "x")
    snp_dist <- make_snp_dist(
        labels,
        list(
            "a-c" = 5,
            "a-d" = 6,
            "b-c" = 7,
            "b-d" = 8,
            "a-x" = 1,
            "b-x" = 9,
            "c-x" = 10,
            "d-x" = 11,
            "a-b" = 3,
            "c-d" = 4
        )
    )
    clusters <- c(a = 1, b = 1, c = 2, d = 2)

    out <- cluster_inter_distances(clusters, snp_dist)

    expect_equal(rownames(out), c("1", "2"))
    # Cluster 1's nearest other-cluster isolate is c (5); its nearest isolate of any kind is x (1).
    expect_equal(out["1", "min_inter_cluster"], 5)
    expect_equal(out["1", "min_inter_isolate"], 1)
    # Cluster 2's nearest other isolate is a (5); x is farther, so both columns are 5.
    expect_equal(out["2", "min_inter_cluster"], 5)
    expect_equal(out["2", "min_inter_isolate"], 5)
})

# --- intra_cluster_duration_metrics --------------------------------------------

test_that("intra_cluster_duration_metrics measures acquisition timing across patients", {
    seqs <- c("g1", "g2", "g3", "g4")
    seq2pt <- c(g1 = "p1", g2 = "p1", g3 = "p2", g4 = "p3")
    dates <- c(g1 = 0, g2 = 5, g3 = 3, g4 = 10) # earliest per patient: p1 = 0, p2 = 3, p3 = 10

    out <- intra_cluster_duration_metrics(seqs, seq2pt, dates)

    expect_equal(unname(out["cluster_start_date"]), 0)
    expect_equal(unname(out["time_to_first_acquisition"]), 3) # 3 - 0
    expect_equal(unname(out["time_to_last_acquisition"]), 10) # 10 - 0
    expect_equal(unname(out["median_time_to_acquisition"]), 6.5) # median(c(3, 10))
})

test_that("intra_cluster_duration_metrics returns zero gaps for a single-patient cluster", {
    out <- intra_cluster_duration_metrics(
        c("g1", "g2"),
        c(g1 = "p1", g2 = "p1"),
        c(g1 = 4, g2 = 9)
    )

    expect_equal(unname(out["cluster_start_date"]), 4)
    expect_equal(unname(out["time_to_first_acquisition"]), 0)
    expect_equal(unname(out["time_to_last_acquisition"]), 0)
    expect_equal(unname(out["median_time_to_acquisition"]), 0)
})

# --- compute_*_by_cluster (integration over a lookup) --------------------------

test_that("compute_genetic_distances_by_cluster tidies one row per multi-patient cluster", {
    lookup <- data.frame(
        cluster = c(1, 1, 2),
        isolate_id = c("a", "b", "c"),
        patient_id = c("p1", "p2", "p3"), # cluster 2 has a single patient -> excluded
        stringsAsFactors = FALSE
    )
    snp_dist <- make_snp_dist(c("a", "b", "c"), list("a-b" = 2, "a-c" = 4, "b-c" = 6))

    out <- compute_genetic_distances_by_cluster(lookup, snp_dist)

    expect_equal(out$cluster, 1)
    expect_equal(out$max_genetic_distance, 2) # only the a-b pair is within cluster 1
})

test_that("compute_duration_metrics_by_cluster tidies one row per multi-patient cluster", {
    lookup <- data.frame(
        cluster = c(1, 1, 2),
        isolate_id = c("g1", "g3", "g4"),
        patient_id = c("p1", "p2", "p3"),
        date = c(0, 3, 10),
        stringsAsFactors = FALSE
    )

    out <- compute_duration_metrics_by_cluster(lookup)

    expect_equal(out$cluster, 1)
    expect_equal(out$time_to_last_acquisition, 3) # cluster 1 patients acquired at days 0 and 3
})
