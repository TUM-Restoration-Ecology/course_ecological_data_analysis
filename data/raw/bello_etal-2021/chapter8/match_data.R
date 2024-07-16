# Matches community, trait, and phylogney data by species names that are
# provided in the form of column names, row names, and tip labels, respectively.
# Options to specify minimum species nr per plot and to scale (i.e. standardize)
# trait values. Community and trait data should be data.frames, not matrices.
# The phylogeny has to be in the phylo format of package ape.
match_data <- function(com = NULL, traits = NULL, tree = NULL, min_richness = 1, 
                       scale = FALSE) {
  #   if(sum(is.null(com), is.null(traits), is.null(tree)) == 1) {
  #     if(!is.null(com)) {
  #       l <- lapply(com, colnames)
  #     }
  #     if(!is.null(traits)) {
  #       l <- lapply(traits, rownames)
  #     }
  #     if(!is.null(tree)) {
  #       l <- lapply(tree, function(x) {x$tip.label})
  #     }
  #   }
  try(com_id <- colnames(com), silent = TRUE)
  try(traits_id <- rownames(traits), silent = TRUE)
  try(tree_id <- tree$tip.label, silent = TRUE)
  
  if(any(c(is.null(com), is.null(traits), is.null(tree)))) {
    if(is.null(com)) l <- list(traits_id, tree_id)
    if(is.null(traits)) l <- list(com_id, tree_id)
    if(is.null(tree)) l <- list(com_id, traits_id)
  } else {
    l <- list(com_id, traits_id, tree_id)
  }
  int <- Reduce(intersect, l)
  try(tree <- drop.tip(tree, setdiff(tree$tip.label, int)), silent = TRUE)
  try(com <- com[, int, drop = FALSE], silent = TRUE)
  if(is.null(tree)) {
    traits <- traits[int, , drop = FALSE]
  } else {
    traits <- traits[tree$tip.label, , drop = FALSE]
  }
  if(min_richness > 1) {
    keep <- which(rowSums(ifelse(com > 0, 1, 0)) > min_richness)
    if(length(keep) == 0) {
      stop(paste("No samples with richness >", min_richness))
    }
    com <- com[keep, ]
    if(any(colSums(com) == 0)) {
      omit <- which(colSums(com) == 0)
      omit_names <- colnames(com)[omit]
      com <- com[, -omit]
      traits <- traits[-omit, ]
      tree <- drop.tip(tree, omit_names)
    }    
  }
  if(scale) {
    names <- rownames(traits)
    traits <- as.data.frame(scale(traits))
    rownames(traits) <- names
  }
  result <- list(com = com, traits = traits, tree = tree)
  result
}

# substitutes the mpd function in the picante packages with a "corrected" version
mpd <- function(samp, dis, abundance.weighted = FALSE) 
{
  N <- dim(samp)[1]
  mpd <- numeric(N)
  for (i in 1:N) {
    sppInSample <- names(samp[i, samp[i, ] > 0])
    if (length(sppInSample) > 1) {
      sample.dis <- dis[sppInSample, sppInSample]
      if (abundance.weighted) {
        sample.weights <- t(as.matrix(samp[i, sppInSample, drop = FALSE])) %*% as.matrix(samp[i, sppInSample, drop = FALSE])
        mpd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis)], sample.weights[lower.tri(sample.weights)])
      }
      else {
        mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
      }
    }
    else {
      mpd[i] <- NA
    }
  }
  mpd
}

