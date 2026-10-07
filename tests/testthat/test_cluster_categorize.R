# Tests for the patient- and cluster-level categorization with small synthetic lookups.

# --- cluster_patient_categorization --------------------------------------------

test_that("cluster_patient_categorization labels an admission-positive index and a convert", {
    lookup <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        patient_id = c("p1", "p2"),
        date = c(0, 10),
        adm_pos = c(TRUE, FALSE),
        prev_surv = c(0, 2), # p2 had a prior screen at day 2, so its day-10 isolate is not the first
        prev_surv_neg = c(FALSE, TRUE),
        stringsAsFactors = FALSE
    )
    surv_df <- data.frame(
        patient_id = c("p1", "p2", "p2"),
        genome_id = c("g1", "g2neg", "g2"),
        surv_date = c(0, 2, 10),
        result = c(1, 0, 1)
    )

    cats <- cluster_patient_categorization(lookup, surv_df)

    expect_equal(cats[["1"]][["p1"]], "index")
    expect_equal(cats[["1"]][["p2"]], "convert")
})

test_that("cluster_patient_categorization labels a not-admission-positive first-surv index as weak", {
    lookup <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        patient_id = c("p1", "p2"),
        date = c(0, 10),
        adm_pos = c(FALSE, FALSE),
        prev_surv = c(0, 2), # p1's first positive is its first surveillance (prev_surv == date)
        prev_surv_neg = c(FALSE, TRUE),
        stringsAsFactors = FALSE
    )
    surv_df <- data.frame(
        patient_id = c("p1", "p2", "p2"),
        genome_id = c("g1", "g2neg", "g2"),
        surv_date = c(0, 2, 10),
        result = c(1, 0, 1)
    )

    cats <- cluster_patient_categorization(lookup, surv_df)
    expect_equal(cats[["1"]][["p1"]], "weak-index")
})

test_that("cluster_patient_categorization labels a culture-only first positive with a prior screen as ambiguous-convert", {
    lookup <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        patient_id = c("p1", "p2"),
        date = c(0, 50),
        adm_pos = c(TRUE, FALSE),
        prev_surv = c(0, 10),
        prev_surv_neg = c(FALSE, FALSE),
        stringsAsFactors = FALSE
    )
    # p2's earliest positive is a culture (genome_id NA) at day 10, after a negative
    # screen at day 2; its sequenced isolate g2 only appears at day 50. The strain at
    # day 10 is unknown, so p2 is ambiguous, and the prior screen makes it a convert.
    surv_df <- data.frame(
        patient_id = c("p1", "p2", "p2", "p2"),
        genome_id = c("g1", NA, NA, "g2"),
        surv_date = c(0, 2, 10, 50),
        result = c(1, 0, 1, 1),
        stringsAsFactors = FALSE
    )

    cats <- cluster_patient_categorization(lookup, surv_df)
    expect_equal(cats[["1"]][["p2"]], "ambiguous-convert")
})

test_that("cluster_patient_categorization labels a culture-only first positive with no prior screen as ambiguous-adm-pos", {
    lookup <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        patient_id = c("p1", "p2"),
        date = c(0, 50),
        adm_pos = c(TRUE, FALSE),
        prev_surv = c(0, 10),
        prev_surv_neg = c(FALSE, FALSE),
        stringsAsFactors = FALSE
    )
    # p2's earliest positive is a culture (genome_id NA) at day 10 with no earlier
    # surveillance; its sequenced isolate g2 only appears at day 50.
    surv_df <- data.frame(
        patient_id = c("p1", "p2", "p2"),
        genome_id = c("g1", NA, "g2"),
        surv_date = c(0, 10, 50),
        result = c(1, 1, 1),
        stringsAsFactors = FALSE
    )

    cats <- cluster_patient_categorization(lookup, surv_df)
    expect_equal(cats[["1"]][["p2"]], "ambiguous-adm-pos")
})

# --- categorize_cluster_overlap ------------------------------------------------

# Surveillance is only consulted in the weak-index branch; a minimal frame suffices elsewhere.
surv_min <- data.frame(
    patient_id = character(0),
    genome_id = character(0),
    surv_date = numeric(0),
    result = numeric(0)
)

test_that("categorize_cluster_overlap returns patient-to-patient for an explained index cluster", {
    lookup <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        patient_id = c("p1", "p2"),
        date = c(0, 10),
        adm_pos = c(TRUE, FALSE),
        prev_surv = c(0, 2),
        prev_surv_neg = c(FALSE, TRUE),
        stringsAsFactors = FALSE
    )
    cluster_overlap_df <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        overlap = c(NA, TRUE)
    )

    cat <- categorize_cluster_overlap(lookup, cluster_overlap_df, surv_min)
    expect_equal(unname(cat["1"]), "patient-to-patient")
})

test_that("categorize_cluster_overlap returns missing-intermediate when a convert is unexplained", {
    lookup <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        patient_id = c("p1", "p2"),
        date = c(0, 10),
        adm_pos = c(TRUE, FALSE),
        prev_surv = c(0, 2),
        prev_surv_neg = c(FALSE, TRUE),
        stringsAsFactors = FALSE
    )
    cluster_overlap_df <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        overlap = c(NA, FALSE)
    )

    cat <- categorize_cluster_overlap(lookup, cluster_overlap_df, surv_min)
    expect_equal(unname(cat["1"]), "missing-intermediate")
})

test_that("categorize_cluster_overlap flags an all-admission-positive cluster", {
    lookup <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        patient_id = c("p1", "p2"),
        date = c(0, 10),
        adm_pos = c(TRUE, TRUE),
        prev_surv = c(0, 0),
        prev_surv_neg = c(FALSE, FALSE),
        stringsAsFactors = FALSE
    )
    cluster_overlap_df <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("g1", "g2"),
        overlap = c(NA, NA)
    )

    cat <- categorize_cluster_overlap(lookup, cluster_overlap_df, surv_min)
    expect_equal(unname(cat["1"]), "all-admission-positive")
})
