library(ape)
library(rotl)
# ---
build_pruned_family_tree <- function(df, family, n_keep, seed) {
  set.seed(seed)

  # 1) species in this family (Genus species strings)
  spp <- unique(df$Scientific[
    df$BLFamilyLatin == family & !is.na(df$Scientific)
  ])
  if (length(spp) < 2) {
    stop(sprintf("Family '%s' has <2 species.", family))
  }

  # 2) resolve to OTT ids (robust: no context)
  tn <- rotl::tnrs_match_names(spp, do_approximate_matching = TRUE)
  ids <- as.numeric(sub("^ott", "", tn$ott_id))
  ids <- unique(ids[!is.na(ids)])
  if (length(ids) < 2) {
    stop("Too few names matched to OTT ids.")
  }

  # helper: call induced_subtree; if a bad node id is reported, drop it and retry once
  get_subtree <- function(ids_vec) {
    tr <- try(rotl::tol_induced_subtree(ott_ids = ids_vec), silent = TRUE)
    if (!inherits(tr, "try-error")) {
      return(tr)
    }
    msg <- conditionMessage(attr(tr, "condition"))
    bad <- as.numeric(sub(".*node_id 'ott([0-9]+)'.*", "\\1", msg))
    if (is.finite(bad)) {
      ids_new <- setdiff(ids_vec, bad)
      if (length(ids_new) >= 2) {
        return(rotl::tol_induced_subtree(ott_ids = ids_new))
      }
    }
    stop(tr)
  }

  # 3) induced subtree -> phylo
  tr <- get_subtree(ids)
  phy <- if (inherits(tr, "phylo")) {
    tr
  } else if (!is.null(tr$tree) && inherits(tr$tree, "phylo")) {
    tr$tree
  } else if (is.character(tr)) {
    ape::read.tree(text = tr)
  } else {
    stop("Unexpected tree format.")
  }

  # 4) clean labels and prune to ≤ n_keep species present in df
  phy$tip.label <- gsub("_", " ", rotl::strip_ott_ids(phy$tip.label))
  keepable <- intersect(phy$tip.label, spp)
  if (!length(keepable)) {
    stop("No overlap between tree tips and species names.")
  }

  if (length(keepable) > n_keep) {
    keepable <- sample(keepable, n_keep)
  }
  phy <- ape::keep.tip(phy, keepable)
  ape::ladderize(phy)
}

# ---

# Align a phylo tree's tips with a trait in BirdBASE
# Returns: list(phy=<pruned phylo>, trait=<named factor in tip order>)
align_traits <- function(
  phy,
  birdbase,
  trait_col,
  species_col,
  level_order = NULL
) {
  stopifnot(inherits(phy, "phylo"))
  stopifnot(trait_col %in% names(birdbase), species_col %in% names(birdbase))

  # pull species & trait, trim, drop missing traits
  sp <- trimws(as.character(birdbase[[species_col]]))
  traw <- trimws(as.character(birdbase[[trait_col]]))
  ok <- !is.na(sp) & nzchar(sp) & !is.na(traw) & nzchar(traw)
  dat <- data.frame(sci = sp[ok], trait = traw[ok], stringsAsFactors = FALSE)

  # collapse duplicates per species: most frequent trait
  if (any(duplicated(dat$sci))) {
    dat <- do.call(
      rbind,
      lapply(split(dat$trait, dat$sci), function(v) {
        tv <- sort(table(v), decreasing = TRUE)
        data.frame(
          sci = names(attr(tv, "dimnames"))[[1]][1], # placeholder, overwritten below
          trait = names(tv)[1],
          stringsAsFactors = FALSE
        )
      })
    )
    dat$sci <- rownames(dat)
    rownames(dat) <- NULL
  }

  # make labels comparable (trees here already use spaces)
  phy$tip.label <- gsub("_", " ", phy$tip.label)
  dat$sci <- gsub("_", " ", dat$sci)

  # keep only overlap; prune tree
  keep <- intersect(phy$tip.label, dat$sci)
  if (length(keep) < 2) {
    stop("Too few overlapping species between tree and BirdBASE.")
  }
  phy <- ape::keep.tip(phy, keep)

  # trait vector in tip order (factor)
  trait_vec <- dat$trait[match(phy$tip.label, dat$sci)]
  trait_fac <- if (is.null(level_order)) {
    factor(trait_vec)
  } else {
    factor(trait_vec, levels = level_order)
  }
  names(trait_fac) <- phy$tip.label

  list(phy = phy, trait = trait_fac)
}
#---
prep_tree_trait <- function(phy, trait, bifurcate = TRUE) {
  stopifnot(inherits(phy, "phylo"), !is.null(names(trait)))

  keep <- intersect(phy$tip.label, names(trait)[!is.na(trait)])
  if (length(keep) < 2L) {
    stop("Too few overlapping tips between tree and trait.")
  }
  phy <- ape::keep.tip(phy, keep)
  tr <- trait[phy$tip.label]
  tr <- droplevels(as.factor(tr))

  if (isTRUE(bifurcate)) {
    phy <- ape::multi2di(phy, random = TRUE)
    phy <- ape::collapse.singles(phy)
  }

  if (is.null(phy$edge.length)) {
    phy <- ape::compute.brlen(phy, method = "Grafen")
  }
  bad <- is.na(phy$edge.length) | phy$edge.length <= 0
  if (any(bad)) {
    phy$edge.length[bad] <- 1e-6
  }

  list(phy = phy, trait = tr)
}


# ---
# kNN sparsifier used in network building
knn_sparsify <- function(M, k, symmetrize = TRUE) {
  n <- nrow(M)
  Sk <- matrix(0, n, n, dimnames = dimnames(M))
  for (i in seq_len(n)) {
    xi <- M[i, ]
    if (all(xi <= 0)) {
      next
    }
    ord <- order(xi, decreasing = TRUE)
    ord <- ord[xi[ord] > 0]
    if (length(ord)) {
      keep <- ord[seq_len(min(k, length(ord)))]
      Sk[i, keep] <- xi[keep]
    }
  }
  if (symmetrize) {
    Sk <- pmax(Sk, t(Sk))
  }
  diag(Sk) <- 0
  Sk
}

# ---
# Core builder: combine phylogeny + diet, sparsify, return nodes/edges
build_edges_nodes <- function(phy, diet, alpha, k) {
  # align diet to tree tips
  diet <- diet[phy$tip.label]

  # phylogenetic similarity from patristic distance
  D <- cophenetic.phylo(phy)
  s0 <- median(D[upper.tri(D)], na.rm = TRUE)
  Sphy <- exp(-D / s0)

  # diet similarity (same class = 1 else 0)
  Sdiet <- outer(diet, diet, FUN = "==") * 1.0
  dimnames(Sdiet) <- list(phy$tip.label, phy$tip.label)

  # combine & sparsify
  S <- alpha * Sphy + (1 - alpha) * Sdiet
  diag(S) <- 0
  S <- knn_sparsify(S, k = k)

  # edge list (collapse to undirected)
  idx <- which(S > 0, arr.ind = TRUE)
  edges <- data.frame(
    source = rownames(S)[idx[, 1]],
    target = colnames(S)[idx[, 2]],
    weight = as.numeric(S[idx]),
    interaction = "phylo_diet",
    stringsAsFactors = FALSE
  )
  edges$pair <- ifelse(
    edges$source < edges$target,
    paste(edges$source, edges$target, sep = "||"),
    paste(edges$target, edges$source, sep = "||")
  )
  edges <- aggregate(weight ~ pair + interaction, data = edges, FUN = max)
  tmp <- do.call(rbind, strsplit(edges$pair, "\\|\\|"))
  edges <- data.frame(
    source = tmp[, 1],
    target = tmp[, 2],
    weight = edges$weight,
    interaction = edges$interaction,
    stringsAsFactors = FALSE
  )

  nodes <- data.frame(
    id = phy$tip.label,
    diet = as.character(diet),
    stringsAsFactors = FALSE
  )

  list(edges = edges, nodes = nodes)
}

# ---
# Send to Cytoscape
send_to_cyto <- function(nodes, edges, title, collection = "Bird_Networks") {
  cytoscapePing()
  createNetworkFromDataFrames(
    nodes = nodes,
    edges = edges,
    title = title,
    collection = collection
  )
}
# ---
