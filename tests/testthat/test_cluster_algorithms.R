# Tests for the two clustering algorithms with small, hand-built toy inputs.
#
# get_tn_clusters_snp_thresh is checked on a SNP distance matrix with two well-separated groups;
# get_tn_clusters_sv_index on a tiny alignment with an outgroup, a clade that shares defining
# variants, and an admission-positive index preceding two converts.

# --- get_tn_clusters_snp_thresh ------------------------------------------------

# Two tight pairs {a, b} and {c, d}, 1 SNP apart within a pair and 20 SNPs apart between pairs.
two_group_snp_dist <- function() {
    labs <- c("a", "b", "c", "d")
    m <- matrix(20, nrow = 4, ncol = 4, dimnames = list(labs, labs))
    diag(m) <- 0
    m["a", "b"] <- m["b", "a"] <- 1
    m["c", "d"] <- m["d", "c"] <- 1
    m
}

test_that("get_tn_clusters_snp_thresh recovers two well-separated groups", {
    clusters <- get_tn_clusters_snp_thresh(two_group_snp_dist(), 5)

    # one assignment per isolate, names preserved
    expect_equal(names(clusters), c("a", "b", "c", "d"))
    expect_length(unique(clusters), 2)
    # each tight pair stays together, and the two pairs are kept apart
    expect_equal(unname(clusters["a"]), unname(clusters["b"]))
    expect_equal(unname(clusters["c"]), unname(clusters["d"]))
    expect_false(unname(clusters["a"]) == unname(clusters["c"]))
})

test_that("get_tn_clusters_snp_thresh merges everything above the largest distance", {
    # a threshold past the 20-SNP gap collapses both groups into a single cluster
    clusters <- get_tn_clusters_snp_thresh(two_group_snp_dist(), 25)
    expect_length(unique(clusters), 1)
})

# --- get_tn_clusters_sv_index --------------------------------------------------

test_that("get_tn_clusters_sv_index groups a defining-variant clade, separate from the outgroup", {
    # Outgroup O first, then a clade {i1, i2, i3} sharing two clade-defining variants
    # (positions 1-2, all C against the outgroup's A) plus one private mutation each. i1 is the
    # admission-positive index (day 0); i2 and i3 are converts (days 10 and 20).
    seqs <- c(
        O = "aaaaaaaaaaaa",
        i1 = "ccaaaaaaaaaa",
        i2 = "ccgaaaaaaaaa",
        i3 = "ccataaaaaaaa"
    )
    mat <- do.call(rbind, strsplit(seqs, ""))
    rownames(mat) <- names(seqs)
    dna_aln <- ape::as.DNAbin(mat)

    snp_dist <- get_snp_dist_matrix(dna_aln)
    tree <- get_phylo_tree(dna_aln, snp_dist, "pars")

    seq2pt <- c(O = "p0", i1 = "p1", i2 = "p2", i3 = "p3")
    dates <- c(O = 0, i1 = 0, i2 = 10, i3 = 20)

    clusters <- suppressMessages(get_tn_clusters_sv_index(
        dna_aln,
        snp_dist,
        adm_seqs = "i1",
        adm_pos_pt_seqs = "i1",
        seq2pt,
        dates,
        tree
    ))

    # one assignment per sequence
    expect_equal(length(clusters), nrow(dna_aln))
    # the three clade members land in a single shared cluster ...
    expect_length(unique(clusters[c("i1", "i2", "i3")]), 1)
    # ... that excludes the outgroup
    expect_false(unname(clusters["O"]) == unname(clusters["i1"]))
})
