# Tests for the spatial/temporal overlap helpers with small synthetic trace matrices.
#
# Trace matrices are patient-by-date (patients in rows, day positions in columns); a positive
# value marks where a patient was on a given day. Overlap is counted over the window from the
# donor's previous surveillance date up to the recipient's collection date.

# A lookup with two patients, each contributing one isolate.
make_lookup <- function(prev_surv = c(1, 4), dates = c(3, 7), adm_pos = c(FALSE, FALSE)) {
    data.frame(
        cluster = c(1, 1),
        isolate_id = c("i1", "i2"),
        patient_id = c("p1", "p2"),
        date = dates,
        adm_pos = adm_pos,
        prev_surv = prev_surv,
        prev_surv_neg = c(TRUE, TRUE),
        stringsAsFactors = FALSE
    )
}

# --- isolate_isolate_overlap (concurrent co-location) ---------------------------

test_that("isolate_isolate_overlap counts concurrent co-location within the donor window", {
    trace <- matrix(0, nrow = 2, ncol = 12, dimnames = list(c("p1", "p2"), as.character(1:12)))
    trace["p1", 1:6] <- 1
    trace["p2", 4:9] <- 1

    out <- isolate_isolate_overlap(make_lookup(), trace)

    # Only i1 -> i2 is valid (for i2 -> i1 the donor window starts after the recipient date).
    # Within columns 1:7 both patients are present at positions 4, 5, 6 -> 3 overlap days.
    expect_equal(nrow(out), 1)
    expect_equal(out$iso_donor, "i1")
    expect_equal(out$iso_recipient, "i2")
    expect_equal(out$overlap_days, 3L)
})

test_that("isolate_isolate_overlap returns an empty, well-formed frame when no pair is valid", {
    # Both isolates from the same patient -> every ordered pair is skipped.
    lookup <- make_lookup()
    lookup$patient_id <- c("p1", "p1")
    trace <- matrix(1, nrow = 1, ncol = 12, dimnames = list("p1", as.character(1:12)))

    out <- isolate_isolate_overlap(lookup, trace)

    expect_equal(nrow(out), 0)
    expect_equal(names(out), c("iso_donor", "iso_recipient", "overlap_days"))
})

test_that("isolate_isolate_overlap drops isolates with missing dates or absent patients", {
    lookup <- make_lookup()
    lookup$date[2] <- NA
    trace <- matrix(0, nrow = 2, ncol = 12, dimnames = list(c("p1", "p2"), as.character(1:12)))
    trace["p1", 1:6] <- 1
    trace["p2", 4:9] <- 1

    # With i2 dropped there is no recipient with a valid date -> no pairs.
    expect_equal(nrow(isolate_isolate_overlap(lookup, trace)), 0)
})

# --- isolate_isolate_sequential_overlap (indirect, shared location over time) ----

test_that("sequential overlap counts later visits to a location the donor occupied earlier", {
    trace <- matrix(0, nrow = 2, ncol = 12, dimnames = list(c("p1", "p2"), as.character(1:12)))
    trace["p1", 1:3] <- 5 # donor in room 5 early
    trace["p2", 5:8] <- 5 # recipient in the same room later, no concurrent overlap

    out <- isolate_isolate_sequential_overlap(
        make_lookup(prev_surv = c(1, 5), dates = c(2, 7)),
        trace
    )

    # i1 -> i2: within columns 1:7 the recipient is in room 5 at positions 5, 6, 7, all after the
    # donor first occupied it (position 1) -> 3 sequential-overlap days.
    expect_equal(out$iso_donor, "i1")
    expect_equal(out$iso_recipient, "i2")
    expect_equal(out$overlap_days, 3L)
})

test_that("sequential overlap is zero when the pair has concurrent co-location", {
    trace <- matrix(0, nrow = 2, ncol = 12, dimnames = list(c("p1", "p2"), as.character(1:12)))
    trace["p1", 1:6] <- 1
    trace["p2", 4:9] <- 1 # overlaps p1 concurrently at positions 4-6

    out <- isolate_isolate_sequential_overlap(make_lookup(), trace)
    expect_equal(out$overlap_days[out$iso_donor == "i1" & out$iso_recipient == "i2"], 0L)
})

# --- cluster_isolate_overlap ---------------------------------------------------

test_that("cluster_isolate_overlap flags recipients with overlap and NAs admission-positives", {
    lookup <- make_lookup(adm_pos = c(TRUE, FALSE))
    iso_overlap_df <- data.frame(
        iso_donor = "i1",
        iso_recipient = "i2",
        overlap_days = 3L,
        stringsAsFactors = FALSE
    )

    out <- cluster_isolate_overlap(lookup, iso_overlap_df)

    expect_equal(out$cluster, c(1, 1))
    expect_equal(out$isolate_id, c("i1", "i2"))
    # i1 is admission-positive -> NA; i2 is a recipient with overlap -> TRUE.
    expect_true(is.na(out$overlap[out$isolate_id == "i1"]))
    expect_true(out$overlap[out$isolate_id == "i2"])
})

# --- fraction_convert_events_with_overlap --------------------------------------

test_that("fraction_convert_events_with_overlap pools convert events and carries counts", {
    lookup <- make_lookup(adm_pos = c(TRUE, FALSE))
    cluster_overlap_df <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("i1", "i2"),
        overlap = c(NA, TRUE)
    )

    frac <- fraction_convert_events_with_overlap(cluster_overlap_df, lookup)

    # p1 is admission-positive (not a convert); p2 is a previous-negative convert with overlap.
    expect_equal(unname(frac["1"]), 1)
    expect_equal(unname(attr(frac, "n_converts")["1"]), 1)
    expect_equal(unname(attr(frac, "n_overlap")["1"]), 1)
})

test_that("fraction_convert_events_with_overlap returns NA for a cluster with no convert events", {
    # Both isolates admission-positive -> no convert events -> NA fraction, zero denominator.
    lookup <- make_lookup(adm_pos = c(TRUE, TRUE))
    lookup$prev_surv_neg <- c(FALSE, FALSE)
    cluster_overlap_df <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("i1", "i2"),
        overlap = c(NA, NA)
    )

    frac <- fraction_convert_events_with_overlap(cluster_overlap_df, lookup)
    expect_true(is.na(frac["1"]))
    expect_equal(unname(attr(frac, "n_converts")["1"]), 0)
})
