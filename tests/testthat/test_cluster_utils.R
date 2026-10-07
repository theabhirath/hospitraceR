# Tests for the pure cluster-bookkeeping helpers with small synthetic examples

test_that("remove_singleton_clusters drops only single-sequence clusters", {
    clusters <- c(a = 1, b = 1, c = 2, d = 3, e = 3)
    out <- remove_singleton_clusters(clusters)

    # cluster 2 is the only singleton and should be removed; names are preserved
    expect_equal(out, c(a = 1, b = 1, d = 3, e = 3))
})

test_that("remove_singleton_clusters returns empty when all clusters are singletons", {
    clusters <- c(a = 1, b = 2, c = 3)
    expect_length(remove_singleton_clusters(clusters), 0)
})

test_that("remove_singleton_clusters handles non-contiguous cluster labels", {
    # cluster 5 is the singleton; its label does not equal its table position
    clusters <- c(a = 1, b = 1, c = 5, d = 8, e = 8)
    out <- remove_singleton_clusters(clusters)

    expect_equal(out, c(a = 1, b = 1, d = 8, e = 8))
})

test_that("get_non_single_patient_clusters keeps clusters with >1 distinct patient", {
    lookup <- data.frame(
        cluster = c(1, 1, 2, 2, 3),
        patient_id = c("p1", "p2", "p3", "p3", "p4"),
        stringsAsFactors = FALSE
    )

    # cluster 1 has two patients; cluster 2 has one patient (sampled twice);
    # cluster 3 is a single isolate -> only cluster 1 qualifies
    expect_equal(get_non_single_patient_clusters(lookup), 1)
})

test_that("flatten_cluster_patient_categorization yields one row per cluster-patient pair", {
    categorization <- list(
        "1" = c(p1 = "index", p2 = "convert"),
        "2" = c(p3 = "index")
    )

    df <- flatten_cluster_patient_categorization(categorization)

    expect_equal(names(df), c("cluster", "patient_id", "category"))
    expect_equal(nrow(df), 3)
    expect_equal(df$cluster, c("1", "1", "2"))
    expect_equal(df$patient_id, c("p1", "p2", "p3"))
    expect_equal(df$category, c("index", "convert", "index"))
})

test_that("flatten_cluster_overlap_categorization yields one row per cluster", {
    categorization <- c("1" = "patient-to-patient", "2" = "inexplicable")

    df <- flatten_cluster_overlap_categorization(categorization)

    expect_equal(names(df), c("cluster", "category"))
    expect_equal(df$cluster, c("1", "2"))
    expect_equal(df$category, c("patient-to-patient", "inexplicable"))
})

# --- get_isolate_lookup --------------------------------------------------------

# A three-sequence DNAbin alignment; the bases are irrelevant to the lookup, only the labels are.
make_dna_aln <- function(labels) {
    mat <- matrix(rep(c("a", "c", "g", "t"), length(labels)), nrow = length(labels), byrow = TRUE)
    rownames(mat) <- labels
    ape::as.DNAbin(mat)
}

test_that("get_isolate_lookup derives prev_surv and prev_surv_neg from surveillance history", {
    dna_aln <- make_dna_aln(c("s1", "s2", "s3"))
    seq2pt <- c(s1 = "p1", s2 = "p2", s3 = "p3")
    dates <- c(s1 = 5, s2 = 20, s3 = 30)
    adm_seqs <- c("s1")
    clusters <- c(s1 = 1, s2 = 1, s3 = 1)
    surv_df <- data.frame(
        patient_id = c("p2", "p3"),
        genome_id = c("x2", "x3"),
        surv_date = c(10, 15), # both strictly before their patient's isolate date
        result = c(0, 1) # p2's prior screen negative; p3's prior screen positive
    )

    lookup <- get_isolate_lookup(clusters, dna_aln, seq2pt, adm_seqs, dates, surv_df)

    expect_equal(
        names(lookup),
        c("isolate_id", "patient_id", "date", "cluster", "adm_pos", "prev_surv", "prev_surv_neg")
    )
    expect_equal(lookup$adm_pos, c(TRUE, FALSE, FALSE))
    # p1 has no surveillance (NA); p2 has a prior screen at day 10; p3 at day 15
    expect_equal(lookup$prev_surv, c(NA, 10, 15))
    # only p2 is a genuine previous-negative -> positive conversion; p3 was already positive
    expect_equal(lookup$prev_surv_neg, c(NA, TRUE, FALSE))
})

test_that("get_isolate_lookup drops isolates missing a patient, date or cluster", {
    dna_aln <- make_dna_aln(c("s1", "s2"))
    seq2pt <- c(s1 = "p1") # s2 has no patient mapping -> dropped by na.omit()
    dates <- c(s1 = 5, s2 = 20)
    clusters <- c(s1 = 1, s2 = 1)
    surv_df <- data.frame(
        patient_id = character(0),
        genome_id = character(0),
        surv_date = numeric(0),
        result = numeric(0)
    )

    lookup <- get_isolate_lookup(clusters, dna_aln, seq2pt, adm_seqs = "s1", dates, surv_df)
    expect_equal(lookup$isolate_id, "s1")
})

# --- remap_cluster_values (internal) -------------------------------------------

test_that("remap_cluster_values numbers values in order of appearance, splitting the special value", {
    out <- hospitraceR:::remap_cluster_values(c(x = 0, y = 5, z = 5, w = 0, v = 3), special_val = 0)

    # 0 is special so each 0 gets its own id; 5 is shared; order of first appearance drives the ids
    expect_equal(out, c(x = 1, y = 2, z = 2, w = 3, v = 4))
})

# --- enforce_monophyly (internal) ----------------------------------------------

test_that("enforce_monophyly expand absorbs intervening isolates into the smallest clade", {
    tree <- ape::read.tree(text = "((a,b),(c,d));")
    clusters <- c(a = 1, b = 1, c = 1, d = 2) # {a,b,c} is not monophyletic on this tree

    out <- suppressMessages(hospitraceR:::enforce_monophyly(clusters, tree, "expand"))

    # the smallest clade containing a, b and c is the whole tree, so d is absorbed: one cluster
    expect_equal(length(unique(out)), 1)
})

test_that("enforce_monophyly break_down splits into the largest pure clades", {
    tree <- ape::read.tree(text = "((a,b),(c,d));")
    clusters <- c(a = 1, b = 1, c = 1, d = 2)

    out <- suppressMessages(hospitraceR:::enforce_monophyly(clusters, tree, "break_down"))

    # {a,b} is a pure clade and stays together; c splits off on its own; d is untouched
    expect_equal(unname(out["a"]), unname(out["b"]))
    expect_false(unname(out["a"]) == unname(out["c"]))
    expect_false(unname(out["c"]) == unname(out["d"]))
    expect_equal(length(unique(out)), 3)
})
