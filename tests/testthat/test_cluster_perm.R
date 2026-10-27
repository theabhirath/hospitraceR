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

    expect_named(res, c("observed", "permuted", "valid_clusters"))

    trace_types <- c("facility", "floor", "room", "seq_facility", "seq_floor", "seq_room")
    # both two-patient clusters are valid
    expect_equal(sort(res$valid_clusters), c(1, 2))

    # observed: one row per trace type, pooled convert-with-overlap and convert counts
    expect_equal(dim(res$observed), c(length(trace_types), 2L))
    expect_equal(rownames(res$observed), trace_types)
    expect_equal(colnames(res$observed), c("n_overlap", "n_converts"))
    # permuted: the same two counts for every permutation
    expect_equal(dim(res$permuted), c(length(trace_types), 2L, nperm))

    # each convert shares facility, floor and room with its own index for 5 days after its
    # negative screen -> 2 of 2 converts explained at every concurrent level
    for (tt in c("facility", "floor", "room")) {
        expect_equal(unname(res$observed[tt, ]), c(2, 2))
    }
    # every convert already has a same-time explanation, so no convert is left for the
    # sequential rows to evaluate
    for (tt in c("seq_facility", "seq_floor", "seq_room")) {
        expect_equal(unname(res$observed[tt, ]), c(0, 0))
    }

    # a permutation can never evaluate more converts than exist
    expect_true(all(res$permuted[, "n_overlap", ] <= res$permuted[, "n_converts", ]))
    expect_true(all(res$permuted[c("facility", "floor", "room"), "n_converts", ] <= 2))
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
