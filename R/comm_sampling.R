#' Function to generate a phylogenetically overdispersed community
#' @param pdist Phylogenetic distance matrix (as generated by \code{\link{cophenetic}})
#' Row names must contain the species names.
#' @param n.spp Number of species to randomly sample from full species pool
#' @param prop Proportion of randomly sampled species to keep at end of overdispersion
#' algorithm
#' @details This function iteratively removed one species from the pair with minimum distance
#' until the number of remaining species equals \code{n.spp*prop}. The initial random sampling
#' makes sure there is some stochasticity to the resulting communities.
#' @export
make_OD_comm <- function(pdist, n.spp, prop) {
  samp <- sample.int(nrow(pdist), n.spp)
  tempdist <- as.matrix(pdist)[samp, samp]
  diag(tempdist) <- max(tempdist) + 1
  while (nrow(tempdist) > ceiling(n.spp*prop)) {
    mins <- apply(tempdist, 1, min)
    themin <- which.min(mins)
    tempdist <- tempdist[-themin, -themin]
  }
  diag(tempdist) <- 0
  return(rownames(tempdist))
}

#' Function to generate a phylogenetically overdispersed communities
#' @param pdist Phylogenetic distance matrix (as generated by \code{\link{cophenetic}}). 
#' Row names must contain the species names.
#' @param n.spp Integer vector specifying the umber of species to randomly sample from 
#' full species pool, for each community.
#' @param prop Numeric vector specifying proportion of randomly sampled species to keep 
#' at end of overdispersion algorithm, for each community.
#' @details This function generates multiple overdispersed communities by iteratively 
#' removed one species from the pair with minimum distance until the number of remaining 
#' species equals \code{n.spp*prop}, where \code{n.spp} and \code{prop} are vectors
#' of length equal to the number of desired communities. The initial random sampling
#' makes sure there is some stochasticity to the resulting communities.
#' @export
make_OD_comms <- function(pdist, n.spp, prop){
  splist <- rownames(pdist)
  comms <- mapply(make_OD_comm, list(pdist = pdist), n.spp, prop, SIMPLIFY = FALSE)
  PAs <- lapply(comms, function(x) (splist %in% x)*1)
  PAmat <- do.call(rbind, PAs)
  colnames(PAmat) <- splist
  rownames(PAmat) <- 1:nrow(PAmat)
  return(PAmat)
}

#' Function to generate a phylogenetically underdispersed community
#' @param pdist Phylogenetic distance matrix (as generated by \code{\link{cophenetic}})
#' Row names must contain the species names.
#' @param n.spp Number of species to randomly sample from full species pool
#' @param prop Proportion of randomly sampled species to keep at end of underdispersion
#' algorithm
#' @details This function starts with a randomly chosen species, and then adds the most 
#' closely related species. It then adds the most closely related species to the newly
#' chosen species, iterating through this procedure until the number of added species equals 
#' \code{n.spp*prop}. The initial random sampling adds extra stochasticity to the resulting 
#' communities.
#' @export
make_UD_comm <- function(pdist, n.spp, prop){
  samp <- sample.int(nrow(pdist), n.spp)
  tempdist <- as.matrix(pdist)[samp, samp]
  diag(tempdist) <- max(tempdist) + 1
  keep <- sample(rownames(tempdist), 1)
  while (length(keep) < ceiling(n.spp*prop)){
    lkeep <- keep[length(keep)]
    themin <- names(which.min(tempdist[lkeep, ]))
    keep <- c(keep, themin)
    rem <- which(rownames(tempdist) %in% lkeep)
    tempdist <- tempdist[-rem, -rem]
  }
  diag(tempdist) <- 0
  return(rownames(pdist[keep, keep]))
}

#' Function to generate a phylogenetically underdispersed communities
#' @param pdist Phylogenetic distance matrix (as generated by \code{\link{cophenetic}}). 
#' Row names must contain the species names.
#' @param n.spp Integer vector specifying the number of species to randomly sample from 
#' full species pool, for each community.
#' @param prop Numeric vector specifying proportion of randomly sampled species to keep 
#' at end of overdispersion algorithm, for each community.
#' @details This function generates multiple underdispersed communities by starting
#' with a randomly chosen species, and then adds the most closely related species. 
#' It then adds the most closely related species to the newly chosen species, iterating 
#' through this procedure until the number of added species equals \code{n.spp*prop}. 
#' The initial random sampling adds extra stochasticity to the resulting communities.
#' @export
make_UD_comms<-function(pdist, n.spp, prop){
  splist <- rownames(pdist)
  comms <- mapply(make_UD_comm, list(pdist = pdist), n.spp, prop, SIMPLIFY = FALSE)
  PAs <- lapply(comms, function(x) (splist %in% x)*1)
  PAmat <- do.call(rbind, PAs)
  colnames(PAmat) <- splist
  rownames(PAmat) <- 1:nrow(PAmat)
  return(PAmat)
}

#' Function to generate a phylogenetically random communities.
#' @param pdist Phylogenetic distance matrix (as generated by \code{\link{cophenetic}}). 
#' Row names must contain the species names.
#' @param n.spp Integer vector specifying the number of species to randomly sample from 
#' full species pool, for each community.
#' @param prop Numeric vector specifying proportion of randomly sampled species to keep 
#' for each community.
#' @details This function generates multiple random communities by  randomly
#' choosing \code{n.spp*prop} species.
#' @export
make_rand_comms <- function(pdist, n.spp, prop) {
  splist <- rownames(pdist)
  comms <- mapply(function(x, y) sample(splist, ceiling(x*y)), n.spp, prop)
  PAs <- lapply(comms, function(x) (splist %in% x)*1)
  PAmat <- do.call(rbind, PAs)
  colnames(PAmat) <- splist
  rownames(PAmat) <- 1:nrow(PAmat)
  return(PAmat)
}
