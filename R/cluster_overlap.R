#' Calculate isolate-isolate co-location overlap
#'
#' @description
#' For every ordered pair of isolates, counts the days on which the two patients were in the same
#' place in `trace_mat`, within the window from the donor's previous surveillance date up to the
#' recipient's collection date. This is the concurrent (spatiotemporal) overlap; for indirect
#' overlap via a shared location at different times, see [isolate_isolate_sequential_overlap()].
#'
#' @param isolate_lookup A lookup table for isolates and their cluster assignments with other
#'                       relevant epidemiological information. See [get_isolate_lookup()].
#' @param trace_mat A patient-by-date matrix of location occupancy (patients in rows, dates in
#'                  columns); a positive value marks a patient's location on a given day.
#'
#' @return A data frame with columns `iso_donor`, `iso_recipient` and `overlap_days`.
#'
#' @export
isolate_isolate_overlap <- function(isolate_lookup, trace_mat) {
    # only keep isolates that:
    # - have non-missing dates and prev_surv
    # - have patients present in the trace matrix
    valid <- !is.na(isolate_lookup$date) &
        !is.na(isolate_lookup$prev_surv) &
        isolate_lookup$patient_id %in% rownames(trace_mat)
    iso <- isolate_lookup[valid, ]
    num_isolates <- nrow(iso)

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
            if (pt_i == pt_j || last_surv >= recipient_date) {
                next
            }
            # get the time range between the donor and recipient
            # days start at 0, trace columns at 1 (named "0", "1", ...)
            cols <- as.integer(last_surv):as.integer(recipient_date) + 1L
            # get the trace matrix subset for the donor and recipient
            pt_sub <- trace_mat[c(row_i, row_j), cols, drop = FALSE]
            # count the number of days where the donor and recipient are both present
            overlap_days <- pt_sub[1, ] == pt_sub[2, ] & pt_sub[1, ] > 0 & pt_sub[2, ] > 0
            # increment the counter and add the donor, recipient, and overlap to the vectors
            k <- k + 1L
            iso_donor_vec[k] <- iso_donor
            iso_rec_vec[k] <- iso_recipient
            overlap_days_vec[k] <- sum(overlap_days)
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
#' shared location. Pairs with any concurrent co-location are excluded, since that is direct
#' overlap and is handled by [isolate_isolate_overlap()].
#'
#' @inheritParams isolate_isolate_overlap
#'
#' @return A data frame with columns `iso_donor`, `iso_recipient` and `overlap_days`.
#'
#' @export
isolate_isolate_sequential_overlap <- function(isolate_lookup, trace_mat) {
    # only keep isolates that:
    # - have non-missing dates and prev_surv
    # - have patients present in the trace matrix
    valid <- !is.na(isolate_lookup$date) &
        !is.na(isolate_lookup$prev_surv) &
        isolate_lookup$patient_id %in% rownames(trace_mat)
    iso <- isolate_lookup[valid, ]
    num_isolates <- nrow(iso)

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
            if (pt_i == pt_j || last_surv > recipient_date) {
                next
            }
            # get the time range between the donor and recipient
            # days start at 0, trace columns at 1 (named "0", "1", ...)
            cols <- as.integer(last_surv):as.integer(recipient_date) + 1L
            # get the trace matrix subset for the donor and recipient
            pt_sub <- trace_mat[c(row_i, row_j), cols, drop = FALSE]
            # if there is ANY spatiotemporal (concurrent) overlap, this pair is
            # invalid for sequential overlap — set to 0
            donor_locs <- pt_sub[1, ]
            recip_locs <- pt_sub[2, ]
            has_concurrent <- any(
                donor_locs == recip_locs & donor_locs > 0 & recip_locs > 0
            )
            seq_overlap <- 0L
            n_days <- length(donor_locs)
            if (!has_concurrent && n_days > 1L) {
                # count days where the recipient is at a location that the donor
                # had previously been at (on a strictly earlier day)
                donor_present <- donor_locs > 0
                unique_donor_locs <- unique(donor_locs[donor_present])
                for (loc in unique_donor_locs) {
                    first_donor_day <- which(donor_locs == loc)[1L]
                    recip_at_loc <- which(recip_locs == loc)
                    seq_overlap <- seq_overlap + sum(recip_at_loc > first_donor_day)
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

#' Summarize isolate overlap within each cluster
#'
#' @description
#' For each isolate in a cluster (treated as a recipient), reports whether any other isolate in the
#' same cluster (as a donor) overlaps with it, using the isolate-pair overlaps from
#' [isolate_isolate_overlap()] or [isolate_isolate_sequential_overlap()]. Admission-positive
#' isolates are reported as `NA`, since they were not acquired in the facility.
#'
#' @param isolate_lookup A lookup table for isolates and their cluster assignments with other
#'                       relevant epidemiological information. See [get_isolate_lookup()].
#' @param iso_overlap_df A data frame of isolate-pair overlaps, from [isolate_isolate_overlap()].
#'
#' @return A data frame with columns `cluster`, `isolate_id` and `overlap`.
#'
#' @export
cluster_isolate_overlap <- function(isolate_lookup, iso_overlap_df) {
    # get non-single patient clusters
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

        ## TODO: check if there is a way to refine this overlap to be more accurate
        # currently all this cares about is if there is any overlap, not the actual overlap days
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

    # combine the data frames for each cluster into a single data frame
    do.call(rbind, out_list)
}

#' Calculate the fraction of convert events with overlap
#'
#' @description
#' Reduces the cluster-level overlaps to, for each cluster, the fraction of convert events that have
#' an overlap explanation.
#'
#' A convert event is an isolate where (1) the patient has no admission-positive
#' isolates, and (2) the patient had a prior negative surveillance before this
#' positive isolate. Criterion (2) is taken from the `prev_surv_neg` flag in the
#' lookup (see [get_isolate_lookup()]): TRUE only when a surveillance exists before
#' the isolate date and none of those prior surveillances were positive. Checking
#' the previous surveillance date alone is not enough, since that prior surveillance
#' may itself have been positive (meaning the patient was already colonized rather
#' than converting). Each qualifying isolate counts as one convert event, and a
#' convert event with an overlap explanation is a convert event with overlap. The
#' per-cluster fraction is `convert_events_with_overlap / convert_events`.
#'
#' The per-cluster counts that make up each fraction are also returned, as the
#' attributes `n_overlap` (convert events with overlap, the numerator) and
#' `n_converts` (convert events, the denominator). These let callers aggregate a
#' pooled, convert-weighted fraction across clusters
#' (`sum(n_overlap) / sum(n_converts)`) instead of averaging the per-cluster
#' fractions, which would ignore how many convert events each cluster contributes.
#'
#' @param cluster_overlap_df A data frame of per-isolate cluster overlap, from
#'                            [cluster_isolate_overlap()].
#' @param isolate_lookup A lookup table for isolates and their cluster assignments with other
#'                        relevant epidemiological information. See [get_isolate_lookup()].
#'
#' @return A named numeric vector with the fraction of convert events with overlap for each
#'   cluster (`NA` for clusters with no convert events), carrying the per-cluster `n_overlap`
#'   and `n_converts` counts as attributes.
#'
#' @importFrom stats setNames
#' @export
fraction_convert_events_with_overlap <- function(cluster_overlap_df, isolate_lookup) {
    clusters <- unique(cluster_overlap_df$cluster)

    # Patients with any admission-positive isolate. A convert event requires the
    # patient to have none of these (they were not colonized at admission).
    adm_pos_patients <- unique(isolate_lookup$patient_id[isolate_lookup$adm_pos])

    # For each cluster, count its convert events (the denominator) and the convert
    # events that have overlap (the numerator). A convert event is an isolate where
    # (1) the patient has no admission-positive isolates, and (2) the patient had a
    # prior negative surveillance before this (positive) isolate -- captured by the
    # prev_surv_neg flag, which is TRUE only when a prior surveillance exists and
    # none of the surveillances before the isolate were positive (a prior positive
    # would mean the patient was already colonized, not a conversion). An event
    # counts as "with overlap" when that isolate has an overlap explanation. Both
    # counts are kept so callers can pool them rather than averaging per-cluster
    # fractions.
    counts <- vapply(
        clusters,
        function(cl) {
            lookup_sub <- isolate_lookup[isolate_lookup$cluster == cl, ]
            convert_events <- !(lookup_sub$patient_id %in% adm_pos_patients) &
                !is.na(lookup_sub$prev_surv_neg) &
                lookup_sub$prev_surv_neg
            convert_isolates <- lookup_sub$isolate_id[convert_events]
            num_converts <- length(convert_isolates)
            num_converts_with_overlap <- sum(
                cluster_overlap_df$overlap[
                    cluster_overlap_df$isolate_id %in% convert_isolates
                ],
                na.rm = TRUE
            )
            c(n_overlap = num_converts_with_overlap, n_converts = num_converts)
        },
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
