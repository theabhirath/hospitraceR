# Integration tests for the overlap permutation test on a small toy dataset. The permutation is
# random, so we check the structure and invariants of the result rather than exact null values.

# Two clusters, each an admission-positive index (day 0) and a previous-negative convert (day 10).
# i1/i2 -> cluster 1 (patients p1, p2); i3/i4 -> cluster 2 (patients p3, p4). Each convert is
# co-located with its own index mid-stay, so the observed overlap fraction is 1.
toy_perm_inputs <- function() {
    isolates <- c("i1", "i2", "i3", "i4")
    dna_aln <- ape::as.DNAbin(
        matrix("a", nrow = 4, ncol = 5, dimnames = list(isolates, NULL))
    )
    days <- as.character(1:12)
    # patient-by-day trace; each index (p1/p3) and its convert (p2/p4) overlap on days 4-8
    mk_trace <- function(codes) {
        m <- matrix(0, nrow = 4, ncol = 12, dimnames = list(c("p1", "p2", "p3", "p4"), days))
        m["p1", 1:8] <- codes[1]
        m["p2", 4:11] <- codes[2]
        m["p3", 1:8] <- codes[3]
        m["p4", 4:11] <- codes[4]
        m
    }
    list(
        clusters = c(i1 = 1, i2 = 1, i3 = 2, i4 = 2),
        dna_aln = dna_aln,
        seq2pt = c(i1 = "p1", i2 = "p2", i3 = "p3", i4 = "p4"),
        adm_seqs = c("i1", "i3"),
        adm_pos_pt_seqs = c("i1", "i3"),
        dates = c(i1 = 0, i2 = 10, i3 = 0, i4 = 10),
        surv_df = data.frame(
            patient_id = c("p1", "p2", "p2", "p3", "p4", "p4"),
            genome_id = c("i1", "p2neg", "i2", "i3", "p4neg", "i4"),
            surv_date = c(0, 2, 10, 0, 2, 10),
            result = c(1, 0, 1, 1, 0, 1),
            stringsAsFactors = FALSE
        ),
        facility_trace = mk_trace(c(1, 1, 1, 1)),
        floor_trace = mk_trace(c(1, 1, 2, 2)),
        room_trace = mk_trace(c(1, 1, 2, 2))
    )
}

test_that("cluster_overlap_perm_test returns a well-formed result", {
    d <- toy_perm_inputs()
    set.seed(42)
    nperm <- 10
    res <- cluster_overlap_perm_test(
        d$clusters,
        d$dna_aln,
        d$seq2pt,
        d$adm_seqs,
        d$adm_pos_pt_seqs,
        d$dates,
        d$surv_df,
        d$facility_trace,
        d$floor_trace,
        d$room_trace,
        nperm = nperm,
        num_cores = 1
    )

    expect_named(
        res,
        c(
            "observed",
            "observed_n_overlap",
            "observed_n_converts",
            "permuted",
            "permuted_n_overlap",
            "permuted_n_converts",
            "valid_clusters"
        )
    )

    trace_types <- c("facility", "floor", "room", "seq_facility", "seq_floor", "seq_room")
    expect_named(res$observed, trace_types)
    # both two-patient clusters are valid
    expect_equal(sort(res$valid_clusters), c(1, 2))

    # the permuted arrays are (n_valid_clusters x 6 trace types x nperm)
    expected_dim <- c(length(res$valid_clusters), length(trace_types), nperm)
    expect_equal(dim(res$permuted), expected_dim)
    expect_equal(dim(res$permuted_n_overlap), expected_dim)
    expect_equal(dim(res$permuted_n_converts), expected_dim)

    # observed fractions are proportions (or NA for clusters with no convert events)
    obs_facility <- res$observed$facility
    expect_true(all(obs_facility >= 0 & obs_facility <= 1, na.rm = TRUE))
    # each convert is co-located with its index at the facility level -> fraction 1 in both clusters
    expect_equal(unname(obs_facility[c("1", "2")]), c(1, 1))

    # each observed fraction is backed by its convert-with-overlap / convert counts
    expect_named(res$observed_n_overlap, trace_types)
    expect_named(res$observed_n_converts, trace_types)
})

test_that("cluster_overlap_perm_test errors when no cluster has more than one patient", {
    d <- toy_perm_inputs()
    # put every isolate in its own cluster -> no multi-patient cluster to test
    singleton_clusters <- setNames(seq_along(d$clusters), names(d$clusters))

    expect_error(
        cluster_overlap_perm_test(
            singleton_clusters,
            d$dna_aln,
            d$seq2pt,
            d$adm_seqs,
            d$adm_pos_pt_seqs,
            d$dates,
            d$surv_df,
            d$facility_trace,
            d$floor_trace,
            d$room_trace,
            nperm = 5,
            num_cores = 1
        ),
        "more than one patient"
    )
})
