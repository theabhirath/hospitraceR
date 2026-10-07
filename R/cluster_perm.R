#' Permutation test for cluster overlap fractions
#'
#' @description
#' Performs a permutation test to assess whether the observed fraction of converts
#' with overlap is significantly different from what would be expected by chance.
#' Tests overlap at facility, floor, and room levels.
#'
#' @param clusters A named numeric vector of cluster assignments.
#' @param dna_aln A DNA alignment object (used for isolate IDs).
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param adm_seqs A vector of sequence IDs which correspond to admission positive sequences.
#' @param adm_pos_pt_seqs A vector of sequence IDs for patients who are admission positive.
#' @param dates A named vector mapping sequence IDs to dates.
#' @param surv_df A data frame with surveillance data.
#' @param facility_trace A matrix with patient IDs as row names and dates as column names.
#' @param floor_trace A matrix with floor-level trace data (same structure as `facility_trace`).
#' @param room_trace A matrix with room-level trace data (same structure as `facility_trace`).
#' @param nperm Number of permutations to perform.
#' @param num_cores Number of cores for parallel processing.
#' @param include_singleton_acq Whether to evaluate all acquisitions, including those not
#'                              classified to be a part of multi-patient clusters.
#'
#' @returns A list containing:
#'     \itemize{
#'     \item `observed`: A matrix with one row per trace type (facility, floor, room,
#'       seq_facility, seq_floor, seq_room) and columns `n_overlap` (converts with overlap) and
#'       `n_converts` (converts), each summed over the tested clusters
#'     \item `permuted`: A numeric array of dimensions (trace type, count, nperm) holding the same
#'       two sums for every permutation
#'     \item `valid_clusters`: A numeric vector of the cluster IDs tested — those with more than one
#'       patient, plus the single-patient acquisition clusters when `include_singleton_acq` is TRUE
#'     }
#'
#'   The pooled, convert-weighted overlap fraction is `n_overlap / n_converts`; summing the counts
#'   across sequence types before dividing pools them further.
#'
#'   A permutation whose cluster assignment strands (capacity too concentrated to place
#'   every patient) stops the whole test with an error.
#'
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmclapply
#' @export
cluster_overlap_perm_test <- function(
    clusters,
    dna_aln,
    seq2pt,
    adm_seqs,
    adm_pos_pt_seqs,
    dates,
    surv_df,
    facility_trace,
    floor_trace,
    room_trace,
    nperm = 1000,
    num_cores = detectCores() - 1,
    include_singleton_acq = FALSE
) {
    # Create observed isolate_lookup
    observed_lookup <- get_isolate_lookup(
        clusters = clusters,
        dna_aln = dna_aln,
        seq2pt = seq2pt,
        adm_seqs = adm_seqs,
        dates = dates,
        surv_df = surv_df
    )

    # Clusters to analyze: multi-patient, plus singleton acquisitions when asked for
    valid_clusters <- overlap_clusters(observed_lookup, include_singleton_acq)

    if (length(valid_clusters) == 0) {
        # Classed so callers can skip an ST with nothing to test without also
        # swallowing real failures.
        stop(errorCondition(
            "No clusters with more than one patient found.",
            class = "no_valid_clusters"
        ))
    }

    # Isolate-isolate overlaps depend only on patient locations over time, not on cluster
    # assignments, so this O(n^2) step runs once and is reused for every permutation. The list
    # names are the trace types of the result.
    traces <- list(facility = facility_trace, floor = floor_trace, room = room_trace)
    iso_overlaps <- c(
        lapply(traces, isolate_isolate_overlap, isolate_lookup = observed_lookup),
        setNames(
            lapply(traces, isolate_isolate_sequential_overlap, isolate_lookup = observed_lookup),
            paste0("seq_", names(traces))
        )
    )

    observed <- pooled_overlap_counts(
        observed_lookup,
        valid_clusters,
        iso_overlaps,
        include_singleton_acq
    )

    draws <- draw_permuted_clusters(clusters, seq2pt, adm_seqs, adm_pos_pt_seqs, nperm, num_cores)
    perm_results <- pbmclapply_or_stop(
        draws,
        function(perm_clust) {
            perm_lookup <- observed_lookup
            perm_lookup$cluster <- perm_clust[perm_lookup$isolate_id]
            pooled_overlap_counts(
                perm_lookup,
                overlap_clusters(perm_lookup, include_singleton_acq),
                iso_overlaps,
                include_singleton_acq
            )
        },
        num_cores
    )

    list(
        observed = observed,
        permuted = simplify2array(perm_results),
        valid_clusters = valid_clusters
    )
}

# =============================================================================
# Helper functions for permutation tests
# =============================================================================

#' Converts with overlap and converts, summed over the given clusters, for each precomputed
#' isolate-isolate overlap table
#' @param iso_overlaps Named list of isolate-pair overlap tables, one per trace type.
#' @return A matrix with one row per trace type and columns `n_overlap`, `n_converts`.
#' @noRd
pooled_overlap_counts <- function(
    isolate_lookup,
    valid_clusters,
    iso_overlaps,
    include_singleton_acq = FALSE
) {
    lookup_filtered <- isolate_lookup[isolate_lookup$cluster %in% valid_clusters, ]
    t(vapply(
        iso_overlaps,
        function(iso_overlap_df) {
            fr <- fraction_convert_events_with_overlap(
                cluster_isolate_overlap(lookup_filtered, iso_overlap_df, include_singleton_acq),
                lookup_filtered
            )
            c(n_overlap = sum(attr(fr, "n_overlap")), n_converts = sum(attr(fr, "n_converts")))
        },
        numeric(2)
    ))
}

#' Create eligibility matrices for permutation
#'
#' Splits the sequences into the three strata the permutation shuffles separately: index
#' patients' sequences in the cluster they started, their sequences in other clusters, and
#' converts. One row per sequence; `comb` keys a patient's sequences within one cluster, which
#' always move together.
#' @noRd
create_eligibility_matrices <- function(clusters, seq2pt, adm_seqs, adm_pos_pt_seqs) {
    seqs <- names(clusters)
    # Sequence ids are names: integer ids would subscript clusters/seq2pt by position.
    adm_seqs <- intersect(as.character(adm_seqs), seqs)
    comb <- setNames(paste(seq2pt[seqs], clusters, sep = "-"), seqs)

    start <- seqs[comb %in% comb[adm_seqs]]
    not_start <- intersect(setdiff(as.character(adm_pos_pt_seqs), start), seqs)
    convert <- setdiff(seqs, c(start, not_start))

    lapply(
        list(index_start = start, index_not_start = not_start, convert = convert),
        function(s) cbind(seq = s, patient = seq2pt[s], cluster = clusters[s], comb = comb[s])
    )
}

#' Assign one stratum's patients to random clusters, preserving how many clusters each patient
#' spans and how many of the stratum's patients each cluster holds
#' @noRd
assign_pt_clusters <- function(elig_mat, cluster_ids, rand_clusters) {
    pairs <- unique(elig_mat[, c("patient", "cluster", "comb"), drop = FALSE])
    capacity <- c(table(factor(pairs[, "cluster"], levels = cluster_ids)))

    # Patients spanning the most clusters go first, while capacity is still spread out
    for (pt in names(sort(c(table(pairs[, "patient"])), decreasing = TRUE))) {
        pt_combs <- pairs[pairs[, "patient"] == pt, "comb"]
        avail <- names(capacity)[capacity > 0]
        # Earlier draws used up the clusters this patient needs: stranded.
        if (length(avail) < length(pt_combs)) {
            stop("cluster assignment stranded (capacity too concentrated to place every patient)")
        }
        cluster_assign <- sample(avail, length(pt_combs))
        capacity[cluster_assign] <- capacity[cluster_assign] - 1

        for (i in seq_along(pt_combs)) {
            seqs_to_assign <- elig_mat[elig_mat[, "comb"] == pt_combs[i], "seq"]
            rand_clusters[seqs_to_assign] <- as.numeric(cluster_assign[i])
        }
    }
    rand_clusters
}

#' Draw `nperm` permuted cluster assignments
#'
#' The null shared by [cluster_overlap_perm_test()] and the analysis-side hotspot tests: patients
#' are reshuffled among clusters within the strata of `create_eligibility_matrices()`, preserving
#' the number of patients per cluster and of clusters per patient, while a patient's sequences in
#' one cluster stay together. A draw that strands (capacity too concentrated to place every
#' patient) is an error.
#' @return A list of `nperm` named cluster vectors.
#' @noRd
draw_permuted_clusters <- function(clusters, seq2pt, adm_seqs, adm_pos_pt_seqs, nperm, num_cores) {
    elig_mats <- create_eligibility_matrices(clusters, seq2pt, adm_seqs, adm_pos_pt_seqs)
    pbmclapply_or_stop(
        seq_len(nperm),
        function(i) assign_permuted_clusters(clusters, elig_mats),
        num_cores
    )
}

#' pbmclapply that raises a failed worker's error
#'
#' With mc.cores > 1, pbmclapply returns a failed worker as a try-error element instead of
#' raising; re-raise the first so callers never index into one.
#' @noRd
pbmclapply_or_stop <- function(X, FUN, num_cores) {
    res <- pbmclapply(X, FUN, mc.cores = num_cores)
    failed <- vapply(res, inherits, logical(1), "try-error")
    if (any(failed)) {
        stop(conditionMessage(attr(res[[which(failed)[1]]], "condition")))
    }
    res
}

#' Assign permuted clusters
#' @noRd
assign_permuted_clusters <- function(clusters, elig_mats) {
    perm_clust <- setNames(rep(-1, length(clusters)), names(clusters))
    cluster_ids <- sort(unique(clusters))

    # Converts first, then index patients who started clusters, then those who didn't
    for (stage in c("convert", "index_start", "index_not_start")) {
        perm_clust <- assign_pt_clusters(elig_mats[[stage]], cluster_ids, perm_clust)
    }

    perm_clust
}
