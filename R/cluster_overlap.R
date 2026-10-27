#' Calculate isolate-isolate co-location overlap
#'
#' @description
#' For every ordered pair of isolates, counts the days on which the two patients were in the same
#' place in `trace_mat`. This is the concurrent (spatiotemporal) overlap; for indirect overlap via a
#' shared location at different times, see [isolate_isolate_sequential_overlap()].
#'
#' The window runs from the later of the donor's previous surveillance date and the recipient's
#' last known negative screen, up to the recipient's collection date. A recipient whose prior
#' surveillance was all negative (`prev_surv_neg`) is known to have been negative on `prev_surv`,
#' so contact before that screen cannot explain the later positive and is not counted; recipients
#' with no known negative screen are not constrained. Both ends are inclusive (contact on the day
#' of the negative screen may follow the swab).
#'
#' A screen can miss a colonization acquired shortly before it, so the recipient's window opens
#' `grace_days` before its last negative screen rather than on it.
#'
#' A pair co-located on fewer than `min_overlap_days` days is reported with `overlap_days = 0`, so
#' that a same-day bed turnover (one patient leaving and another arriving on the same day, which
#' the daily trace records as both being there) is not taken as contact.
#'
#' @param isolate_lookup A lookup table for isolates and their cluster assignments with other
#'                       relevant epidemiological information. See [get_isolate_lookup()].
#' @param trace_mat A patient-by-date matrix of location occupancy (patients in rows, dates in
#'                  columns); a positive value marks a patient's location on a given day.
#' @param min_overlap_days Minimum number of co-located days in the window for a pair to count as
#'                         overlapping.
#' @param grace_days Days before the recipient's last negative screen that still count, allowing
#'                   for a recent colonization the screen did not yet detect.
#'
#' @return A data frame with columns `iso_donor`, `iso_recipient` and `overlap_days`.
#'
#' @export
isolate_isolate_overlap <- function(
    isolate_lookup,
    trace_mat,
    min_overlap_days = 2L,
    grace_days = 3L
) {
    # only keep isolates that:
    # - have non-missing dates and prev_surv
    # - have patients present in the trace matrix
    valid <- !is.na(isolate_lookup$date) &
        !is.na(isolate_lookup$prev_surv) &
        isolate_lookup$patient_id %in% rownames(trace_mat)
    iso <- isolate_lookup[valid, ]
    num_isolates <- nrow(iso)
    rec_start <- recipient_window_start(iso, grace_days)

    # precompute the patient row indices in trace_mat (avoid repeated character lookups)
    pt_row_idx <- match(as.character(iso$patient_id), rownames(trace_mat))

    # preallocate the maximum possible size (num_isolates^2), then trim at the end
    max_pairs <- num_isolates * num_isolates
    iso_donor_vec <- character(max_pairs)
    iso_rec_vec <- character(max_pairs)
    overlap_days_vec <- integer(max_pairs)
    k <- 0L

    # iterate over all pairs of isolates
    for (i in seq_len(num_isolates)) {
        pt_i <- iso$patient_id[i]
        last_surv <- iso$prev_surv[i] # donor window start
        row_i <- pt_row_idx[i]
        iso_donor <- iso$isolate_id[i]

        for (j in seq_len(num_isolates)) {
            iso_recipient <- iso$isolate_id[j]
            pt_j <- iso$patient_id[j]
            recipient_date <- iso$date[j] # recipient collection date
            row_j <- pt_row_idx[j]
            # window start: no earlier than the recipient's last known negative screen, less the
            # grace period
            window_start <- max(last_surv, rec_start[j])
            # skip if:
            # - the donor and recipient are the same patient
            # - the window starts on or after the recipient collection date
            if (pt_i == pt_j || window_start >= recipient_date) {
                next
            }
            # get the time range between the donor and recipient
            # days start at 0, trace columns at 1 (named "0", "1", ...)
            cols <- as.integer(window_start):as.integer(recipient_date) + 1L
            # get the trace matrix subset for the donor and recipient
            pt_sub <- trace_mat[c(row_i, row_j), cols, drop = FALSE]
            # count the number of days where the donor and recipient are both present
            overlap_days <- sum(pt_sub[1, ] == pt_sub[2, ] & pt_sub[1, ] > 0 & pt_sub[2, ] > 0)
            # increment the counter and add the donor, recipient, and overlap to the vectors
            k <- k + 1L
            iso_donor_vec[k] <- iso_donor
            iso_rec_vec[k] <- iso_recipient
            overlap_days_vec[k] <- if (overlap_days >= min_overlap_days) overlap_days else 0L
        }
    }

    # if no overlaps were found, return an empty data frame
    # otherwise, return the overlaps found as a data frame
    overlap_df <- if (k == 0L) {
        data.frame(
            iso_donor = character(0),
            iso_recipient = character(0),
            overlap_days = integer(0)
        )
    } else {
        data.frame(
            iso_donor = iso_donor_vec[seq_len(k)],
            iso_recipient = iso_rec_vec[seq_len(k)],
            overlap_days = overlap_days_vec[seq_len(k)]
        )
    }

    overlap_df
}

#' Calculate isolate-isolate sequential overlap
#'
#' @description
#' For every ordered pair of isolates, counts the days on which the recipient was at a location the
#' donor had occupied on a strictly earlier day, capturing potential indirect transmission via a
#' shared location. The donor's occupancy counts from the donor's previous surveillance date (the
#' location may stay contaminated after the donor leaves), but the recipient's days count only from
#' its last known negative screen (less `grace_days`), as in [isolate_isolate_overlap()].
#'
#' This is a pure pair-level count: pairs that were also concurrently co-located are not zeroed
#' here. Whether a recipient's acquisition already has a same-time explanation is a property of the
#' recipient and its whole cluster, not of one pair, so that exclusion is applied at the cluster
#' level (see [cluster_overlap_perm_test()]).
#'
#' @inheritParams isolate_isolate_overlap
#'
#' @return A data frame with columns `iso_donor`, `iso_recipient` and `overlap_days`.
#'
#' @export
isolate_isolate_sequential_overlap <- function(isolate_lookup, trace_mat, grace_days = 3L) {
    # only keep isolates that:
    # - have non-missing dates and prev_surv
    # - have patients present in the trace matrix
    valid <- !is.na(isolate_lookup$date) &
        !is.na(isolate_lookup$prev_surv) &
        isolate_lookup$patient_id %in% rownames(trace_mat)
    iso <- isolate_lookup[valid, ]
    num_isolates <- nrow(iso)
    rec_start <- recipient_window_start(iso, grace_days)

    # precompute the patient row indices in trace_mat (avoid repeated character lookups)
    pt_row_idx <- match(as.character(iso$patient_id), rownames(trace_mat))

    # preallocate the maximum possible size (num_isolates^2), then trim at the end
    max_pairs <- num_isolates * num_isolates
    iso_donor_vec <- character(max_pairs)
    iso_rec_vec <- character(max_pairs)
    overlap_days_vec <- integer(max_pairs)
    k <- 0L

    # iterate over all pairs of isolates
    for (i in seq_len(num_isolates)) {
        pt_i <- iso$patient_id[i]
        last_surv <- iso$prev_surv[i] # donor window start
        row_i <- pt_row_idx[i]
        iso_donor <- iso$isolate_id[i]

        for (j in seq_len(num_isolates)) {
            iso_recipient <- iso$isolate_id[j]
            pt_j <- iso$patient_id[j]
            recipient_date <- iso$date[j] # recipient collection date
            row_j <- pt_row_idx[j]
            # skip if:
            # - the donor and recipient are the same patient
            # - the donor window starts after the recipient collection date
            # - the recipient's last known negative is on/after its collection date
            if (pt_i == pt_j || last_surv > recipient_date || rec_start[j] >= recipient_date) {
                next
            }
            # get the time range between the donor and recipient
            # days start at 0, trace columns at 1 (named "0", "1", ...)
            cols <- as.integer(last_surv):as.integer(recipient_date) + 1L
            # get the trace matrix subset for the donor and recipient
            pt_sub <- trace_mat[c(row_i, row_j), cols, drop = FALSE]
            donor_locs <- pt_sub[1, ]
            recip_locs <- pt_sub[2, ]
            # recipient days before its last known negative screen (less the grace period) cannot
            # explain the positive
            # (column 1 is day `last_surv`; rec_start is -Inf when there is no known negative)
            recip_from <- max(1, rec_start[j] - last_surv + 1)
            seq_overlap <- 0L
            if (length(donor_locs) > 1L) {
                # count days where the recipient is at a location that the donor
                # had previously been at (on a strictly earlier day)
                for (loc in unique(donor_locs[donor_locs > 0])) {
                    first_donor_day <- which(donor_locs == loc)[1L]
                    recip_at_loc <- which(recip_locs == loc)
                    seq_overlap <- seq_overlap +
                        sum(recip_at_loc > first_donor_day & recip_at_loc >= recip_from)
                }
            }
            # increment the counter and add the donor, recipient, and overlap to the vectors
            k <- k + 1L
            iso_donor_vec[k] <- iso_donor
            iso_rec_vec[k] <- iso_recipient
            overlap_days_vec[k] <- seq_overlap
        }
    }

    # if no overlaps were found, return an empty data frame
    # otherwise, return the overlaps found as a data frame
    overlap_df <- if (k == 0L) {
        data.frame(
            iso_donor = character(0),
            iso_recipient = character(0),
            overlap_days = integer(0)
        )
    } else {
        data.frame(
            iso_donor = iso_donor_vec[seq_len(k)],
            iso_recipient = iso_rec_vec[seq_len(k)],
            overlap_days = overlap_days_vec[seq_len(k)]
        )
    }

    overlap_df
}

#' Earliest day a recipient's exposure can explain its isolate
#'
#' The recipient's last known negative screen less `grace_days`: `prev_surv - grace_days` when
#' every prior surveillance was negative (`prev_surv_neg`), otherwise `-Inf` (no negative result
#' rules any contact out).
#' @noRd
recipient_window_start <- function(iso, grace_days) {
    known_neg <- !is.na(iso$prev_surv_neg) & iso$prev_surv_neg
    ifelse(known_neg, iso$prev_surv - grace_days, -Inf)
}

#' Summarize isolate overlap within each cluster
#'
#' @description
#' For each isolate in a cluster (treated as a recipient), reports whether any other isolate in the
#' same cluster (as a donor) overlaps with it, using the isolate-pair overlaps from
#' [isolate_isolate_overlap()] or [isolate_isolate_sequential_overlap()]. Admission-positive
#' isolates are reported as `NA`, since they were not acquired in the facility.
#'
#' Only multi-patient clusters are reported by default; set `include_singleton_acq` to also account
#' for acquisitions sitting in single-patient clusters.
#'
#' @param isolate_lookup A lookup table for isolates and their cluster assignments with other
#'                       relevant epidemiological information. See [get_isolate_lookup()].
#' @param iso_overlap_df A data frame of isolate-pair overlaps, from [isolate_isolate_overlap()].
#' @param include_singleton_acq Whether to evaluate all acquisitions, including those not
#'                              classified to be a part of multi-patient clusters.
#'
#' @return A data frame with columns `cluster`, `isolate_id` and `overlap`.
#'
#' @export
cluster_isolate_overlap <- function(isolate_lookup, iso_overlap_df, include_singleton_acq = FALSE) {
    # get non-single patient clusters (the singleton acquisitions are appended below)
    clusters <- get_non_single_patient_clusters(isolate_lookup)

    # list of data frames for each cluster
    out_list <- lapply(clusters, function(cl) {
        idx_cl <- isolate_lookup$cluster == cl
        isolates <- isolate_lookup$isolate_id[idx_cl]
        adm_pos <- isolate_lookup$adm_pos[idx_cl]

        # overlaps within this cluster only
        od <- iso_overlap_df[
            iso_overlap_df$iso_donor %in% isolates & iso_overlap_df$iso_recipient %in% isolates,
            ,
            drop = FALSE
        ]

        # check if there is any overlap for each recipient isolate in the cluster
        od_overlap_days <- od$overlap_days
        any_overlap_by_rec <- if (nrow(od) > 0 && max(od_overlap_days) > 0) {
            # this could be min if we want to look at bare minimum overlap
            tapply(od_overlap_days[od_overlap_days > 0], od$iso_recipient[od_overlap_days > 0], any)
        } else {
            FALSE
        }

        # align to isolates; missing => FALSE; admission-positive => NA
        overlaps <- any_overlap_by_rec[isolates]
        overlaps[is.na(overlaps)] <- FALSE
        overlaps[adm_pos] <- NA

        data.frame(
            cluster = rep(cl, length(isolates)),
            isolate_id = isolates,
            overlap = overlaps,
            row.names = NULL
        )
    })

    if (include_singleton_acq) {
        # Every within-cluster pair of a single-patient cluster is a same-patient pair, and
        # isolate_isolate_overlap() never reports those, so the loop above would find no overlap for
        # these isolates however long it scanned. Emit FALSE directly instead
        idx_singleton <- isolate_lookup$cluster %in% get_singleton_acq_clusters(isolate_lookup)
        out_list <- c(
            out_list,
            list(data.frame(
                cluster = isolate_lookup$cluster[idx_singleton],
                isolate_id = isolate_lookup$isolate_id[idx_singleton],
                # rep() rather than a bare FALSE so the no-singletons case stays a 0-row frame
                # instead of tripping data.frame()'s recycling check
                overlap = rep(FALSE, sum(idx_singleton)),
                row.names = NULL
            ))
        )
    }

    # combine the data frames for each cluster into a single data frame
    do.call(rbind, out_list)
}

#' Calculate the fraction of convert events with overlap
#'
#' @description
#' Reduces the cluster-level overlaps to, for each cluster, the fraction of convert events that have
#' an overlap explanation.
#'
#' A convert event is counted **per patient**: one patient's acquisition within one cluster. A
#' convert isolate is one where (1) the patient has no admission-positive isolates, and (2) the
#' patient had a prior negative surveillance before this positive isolate. Criterion (2) is taken
#' from the `prev_surv_neg` flag in the lookup (see [get_isolate_lookup()]): TRUE only when a
#' surveillance exists before the isolate date and none of those prior surveillances were
#' positive. Checking the previous surveillance date alone is not enough, since that prior
#' surveillance may itself have been positive (meaning the patient was already colonized rather
#' than converting). A patient with several convert isolates in a cluster (e.g. two genomes from one
#' culture) is still one convert event, so sequencing depth does not weight the count; the event has
#' an overlap explanation when any of those isolates does. Per cluster, a patient whose isolates
#' fall in two clusters appears in both; for totals across clusters use [pooled_convert_overlap()],
#' which counts each patient once. The per-cluster fraction is
#' `convert_events_with_overlap / convert_events`. Convert isolates whose `overlap` is `NA` in
#' `cluster_overlap_df` were excluded from the evaluation; a patient all of whose convert isolates
#' are `NA` is left out of both counts.
#'
#' The per-cluster counts that make up each fraction are also returned, as the
#' attributes `n_overlap` (convert events with overlap, the numerator) and
#' `n_converts` (convert events, the denominator). To pool across clusters, use
#' [pooled_convert_overlap()] rather than averaging the per-cluster fractions (which
#' would ignore how many convert events each cluster contributes) or summing these
#' counts (which would count a patient spanning two clusters twice).
#'
#' @param cluster_overlap_df A data frame of per-isolate cluster overlap, from
#'                            [cluster_isolate_overlap()].
#' @param isolate_lookup A lookup table for isolates and their cluster assignments with other
#'                        relevant epidemiological information. See [get_isolate_lookup()].
#'
#' @return A named numeric vector with the fraction of convert events (patients) with overlap for
#'   each cluster (`NA` for clusters with no convert events), carrying the per-cluster `n_overlap`
#'   and `n_converts` counts as attributes.
#'
#' @importFrom stats setNames
#' @export
fraction_convert_events_with_overlap <- function(cluster_overlap_df, isolate_lookup) {
    clusters <- unique(cluster_overlap_df$cluster)
    conv <- convert_isolate_overlap(cluster_overlap_df, isolate_lookup)

    # For each cluster, count its convert patients (the denominator) and the convert patients with
    # an overlap explanation (the numerator). Both counts are kept so callers can pool them rather
    # than averaging per-cluster fractions.
    counts <- vapply(
        clusters,
        function(cl) count_convert_patients(conv[conv$cluster == cl, , drop = FALSE]),
        numeric(2)
    )

    cluster_chr <- as.character(clusters)
    n_overlap <- setNames(counts["n_overlap", ], cluster_chr)
    n_converts <- setNames(counts["n_converts", ], cluster_chr)

    # Per-cluster fraction; NA when a cluster has no convert events.
    fractions <- setNames(
        ifelse(n_converts > 0, n_overlap / n_converts, NA_real_),
        cluster_chr
    )

    # Attach the underlying counts so downstream aggregation can pool them.
    attr(fractions, "n_overlap") <- n_overlap
    attr(fractions, "n_converts") <- n_converts
    fractions
}

#' Pooled convert patients with overlap across clusters
#'
#' @description
#' The pooled counterpart of [fraction_convert_events_with_overlap()]: the number of convert
#' patients with an overlap explanation (`n_overlap`) and of convert patients (`n_converts`) over
#' every cluster in `cluster_overlap_df` together. Each patient counts **once** however many
#' clusters their convert isolates fall in (e.g. genomes from one culture split across two
#' clusters), and is explained when any of those isolates has an overlap explanation. Summing the
#' per-cluster counts of [fraction_convert_events_with_overlap()] instead would count such a
#' patient once per cluster. As there, a patient all of whose convert isolates have `NA` overlap
#' (excluded from the evaluation) is left out of both counts.
#'
#' @inheritParams fraction_convert_events_with_overlap
#'
#' @return A named numeric vector with elements `n_overlap` and `n_converts`.
#'
#' @export
pooled_convert_overlap <- function(cluster_overlap_df, isolate_lookup) {
    count_convert_patients(convert_isolate_overlap(cluster_overlap_df, isolate_lookup))
}

#' Overlap of every convert isolate in the clusters of `cluster_overlap_df`
#'
#' A convert isolate is one whose patient has no admission-positive isolates and had a prior
#' negative surveillance (`prev_surv_neg`, TRUE only when a prior surveillance exists and none
#' before the isolate were positive; a prior positive means the patient was already colonized). Its
#' `overlap` is `NA` when the isolate was excluded from the evaluation (e.g. the sequential test
#' skips converts that already have a same-time explanation) and `FALSE` when the isolate is absent
#' from `cluster_overlap_df`.
#' @return A data frame with columns `cluster`, `patient_id` and `overlap`.
#' @noRd
convert_isolate_overlap <- function(cluster_overlap_df, isolate_lookup) {
    adm_pos_patients <- unique(isolate_lookup$patient_id[isolate_lookup$adm_pos])
    is_convert <- isolate_lookup$cluster %in%
        cluster_overlap_df$cluster &
        !(isolate_lookup$patient_id %in% adm_pos_patients) &
        !is.na(isolate_lookup$prev_surv_neg) &
        isolate_lookup$prev_surv_neg
    idx <- match(isolate_lookup$isolate_id[is_convert], cluster_overlap_df$isolate_id)
    overlap <- cluster_overlap_df$overlap[idx]
    overlap[is.na(idx)] <- FALSE
    data.frame(
        cluster = isolate_lookup$cluster[is_convert],
        patient_id = isolate_lookup$patient_id[is_convert],
        overlap = overlap
    )
}

#' Convert patients with overlap / convert patients, from `convert_isolate_overlap()` rows
#'
#' A patient is evaluated if any of their convert isolates was (non-NA overlap) and explained if
#' any of them has overlap.
#' @noRd
count_convert_patients <- function(conv) {
    c(
        n_overlap = length(unique(conv$patient_id[conv$overlap %in% TRUE])),
        n_converts = length(unique(conv$patient_id[!is.na(conv$overlap)]))
    )
}
