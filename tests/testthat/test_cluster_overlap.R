# Tests for the spatial/temporal overlap helpers with small synthetic trace matrices.
#
# Trace matrices are patient-by-date (patients in rows, day positions in columns); a positive
# value marks where a patient was on a given day, and day d is column d + 1. Overlap is counted
# over the window from the later of the donor's previous surveillance date and the recipient's
# last known negative screen, up to the recipient's collection date.

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

test_that("isolate_isolate_overlap counts concurrent co-location within the window", {
    trace <- matrix(0, nrow = 2, ncol = 12, dimnames = list(c("p1", "p2"), as.character(1:12)))
    trace["p1", 1:6] <- 1
    trace["p2", 4:9] <- 1

    out <- isolate_isolate_overlap(make_lookup(), trace, grace_days = 0L)

    # Only i1 -> i2 is valid (for i2 -> i1 the window starts after the recipient date). The window
    # is days 4 (i2's negative screen) to 7, columns 5:8, where both are present in columns 5, 6.
    expect_equal(nrow(out), 1)
    expect_equal(out$iso_donor, "i1")
    expect_equal(out$iso_recipient, "i2")
    expect_equal(out$overlap_days, 2L)
})

test_that("isolate_isolate_overlap opens the window grace_days before the negative screen", {
    trace <- matrix(0, nrow = 2, ncol = 12, dimnames = list(c("p1", "p2"), as.character(1:12)))
    trace["p1", 1:6] <- 1
    trace["p2", 4:9] <- 1

    # Default 3-day grace: the window opens on day 4 - 3 = 1 (no earlier than the donor's day 1),
    # columns 2:8, so the shared columns 4, 5, 6 all count.
    out <- isolate_isolate_overlap(make_lookup(), trace)
    expect_equal(out$overlap_days[out$iso_donor == "i1" & out$iso_recipient == "i2"], 3L)
})

test_that("isolate_isolate_overlap ignores contact before the recipient's last negative screen", {
    trace <- matrix(0, nrow = 2, ncol = 12, dimnames = list(c("p1", "p2"), as.character(1:12)))
    trace["p1", 1:4] <- 1
    trace["p2", 1:9] <- 1 # together on days 0-3, all before i2's negative screen on day 4

    lookup <- make_lookup(prev_surv = c(0, 4))
    out <- isolate_isolate_overlap(lookup, trace, grace_days = 0L)
    expect_equal(out$overlap_days[out$iso_donor == "i1" & out$iso_recipient == "i2"], 0L)

    # The default 3-day grace reopens days 1-3, which count.
    out <- isolate_isolate_overlap(lookup, trace)
    expect_equal(out$overlap_days[out$iso_donor == "i1" & out$iso_recipient == "i2"], 3L)

    # Without a known negative screen the recipient is not constrained: days 0-3 count.
    lookup$prev_surv_neg <- c(TRUE, FALSE)
    out <- isolate_isolate_overlap(lookup, trace)
    expect_equal(out$overlap_days[out$iso_donor == "i1" & out$iso_recipient == "i2"], 4L)
})

test_that("isolate_isolate_overlap needs min_overlap_days co-located days", {
    trace <- matrix(0, nrow = 2, ncol = 12, dimnames = list(c("p1", "p2"), as.character(1:12)))
    trace["p1", 1:6] <- 1
    trace["p2", 6:9] <- 1 # a single shared day (column 6): a same-day bed turnover

    out <- isolate_isolate_overlap(make_lookup(), trace)
    expect_equal(out$overlap_days, 0L)
    out_1day <- isolate_isolate_overlap(make_lookup(), trace, min_overlap_days = 1)
    expect_equal(out_1day$overlap_days, 1L)
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
    trace["p1", 1:3] <- 5 # donor in room 5 on days 0-2
    trace["p2", 5:8] <- 5 # recipient in the same room on days 4-7, no concurrent overlap

    lookup <- make_lookup(prev_surv = c(1, 5), dates = c(2, 7))
    out <- isolate_isolate_sequential_overlap(lookup, trace, grace_days = 0L)

    # i1 -> i2: the donor occupied room 5 from day 1 (its window start). The recipient was there on
    # days 4-7, but day 4 precedes its negative screen on day 5 -> days 5, 6, 7 count.
    expect_equal(out$iso_donor, "i1")
    expect_equal(out$iso_recipient, "i2")
    expect_equal(out$overlap_days, 3L)

    # With the default 3-day grace the recipient's days count from day 2, so day 4 counts too.
    out <- isolate_isolate_sequential_overlap(lookup, trace)
    expect_equal(out$overlap_days, 4L)
})

test_that("sequential overlap is a pair-level count that does not zero concurrent pairs", {
    trace <- matrix(0, nrow = 2, ncol = 12, dimnames = list(c("p1", "p2"), as.character(1:12)))
    trace["p1", 1:6] <- 1
    trace["p2", 4:9] <- 1 # concurrent on days 3-5, and in the donor's location after day 1

    out <- isolate_isolate_sequential_overlap(make_lookup(), trace, grace_days = 0L)
    # window days 1-7; recipient days from its negative screen (day 4) on: days 4-7 -> 4 days
    expect_equal(out$overlap_days[out$iso_donor == "i1" & out$iso_recipient == "i2"], 4L)
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

test_that("fraction_convert_events_with_overlap leaves NA-overlap converts out of both counts", {
    lookup <- make_lookup(adm_pos = c(FALSE, FALSE))
    # i1 was excluded from this evaluation (NA); i2 is evaluated and has no overlap.
    cluster_overlap_df <- data.frame(
        cluster = c(1, 1),
        isolate_id = c("i1", "i2"),
        overlap = c(NA, FALSE)
    )

    frac <- fraction_convert_events_with_overlap(cluster_overlap_df, lookup)
    expect_equal(unname(attr(frac, "n_converts")["1"]), 1)
    expect_equal(unname(attr(frac, "n_overlap")["1"]), 0)
})

test_that("fraction_convert_events_with_overlap counts converts per patient, not per isolate", {
    # p2 contributes two convert isolates (e.g. two genomes from one culture): one explained.
    lookup <- data.frame(
        cluster = c(1, 1, 1, 1),
        isolate_id = c("i1", "i2", "i2b", "i3"),
        patient_id = c("p1", "p2", "p2", "p3"),
        date = c(3, 7, 7, 8),
        adm_pos = c(TRUE, FALSE, FALSE, FALSE),
        prev_surv = c(1, 4, 4, 5),
        prev_surv_neg = c(FALSE, TRUE, TRUE, TRUE),
        stringsAsFactors = FALSE
    )
    cluster_overlap_df <- data.frame(
        cluster = c(1, 1, 1, 1),
        isolate_id = c("i1", "i2", "i2b", "i3"),
        overlap = c(NA, TRUE, FALSE, FALSE)
    )

    frac <- fraction_convert_events_with_overlap(cluster_overlap_df, lookup)
    # Two convert patients (p2, p3); p2 is explained through i2 -> 1 of 2, not 1 of 3 isolates.
    expect_equal(unname(attr(frac, "n_converts")["1"]), 2)
    expect_equal(unname(attr(frac, "n_overlap")["1"]), 1)

    # A patient whose convert isolates were all excluded (NA) leaves the denominator.
    cluster_overlap_df$overlap[cluster_overlap_df$isolate_id %in% c("i2", "i2b")] <- NA
    frac <- fraction_convert_events_with_overlap(cluster_overlap_df, lookup)
    expect_equal(unname(attr(frac, "n_converts")["1"]), 1)
    expect_equal(unname(attr(frac, "n_overlap")["1"]), 0)
})

test_that("pooled_convert_overlap counts a patient spanning two clusters once", {
    # p2's genomes from one culture fall in clusters 1 and 2; only the cluster-2 one is explained.
    lookup <- data.frame(
        cluster = c(1, 1, 2, 2, 1),
        isolate_id = c("i1", "i2", "i2b", "i4", "i3"),
        patient_id = c("p1", "p2", "p2", "p4", "p3"),
        date = c(3, 7, 7, 2, 8),
        adm_pos = c(TRUE, FALSE, FALSE, TRUE, FALSE),
        prev_surv = c(1, 4, 4, 1, 5),
        prev_surv_neg = c(FALSE, TRUE, TRUE, FALSE, TRUE),
        stringsAsFactors = FALSE
    )
    cluster_overlap_df <- data.frame(
        cluster = c(1, 1, 1, 2, 2),
        isolate_id = c("i1", "i2", "i3", "i2b", "i4"),
        overlap = c(NA, FALSE, FALSE, TRUE, NA)
    )

    # Per cluster, p2 appears in both clusters...
    frac <- fraction_convert_events_with_overlap(cluster_overlap_df, lookup)
    expect_equal(sum(attr(frac, "n_converts")), 3)
    # ...but pooled it is one convert patient, explained through i2b: p2 and p3 -> 1 of 2.
    pooled <- pooled_convert_overlap(cluster_overlap_df, lookup)
    expect_equal(pooled, c(n_overlap = 1, n_converts = 2))
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
