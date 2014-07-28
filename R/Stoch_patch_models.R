#' Stochastic Patch model with patch "types"
#' @param t Timestep. Not used in the simulation model, but required for using \code{\link{deSolve}}
#' @param y Current state of the model. This is a vector with length equal to the number of patches in the site.
#' Each element is an integer that indexes the "species" that occupies that patch. 0 represents an empty patch.
#' @param params A list of parameter for the simulation. The parameters are as follows:
#' \itemize{
#'  \item{"d" - The death rate. Rate at which occupied patches become unoccupied}
#'  \item{"nspec" - The number of "species" in the total species pool}
#'  \item{"poolprob" - The probability that an empty patch will be colonized by a species from the pool}
#'  \item{"ntype" - The number of different patch "types"}
#'  \item{"occmat" - A \code{nspec} by \code{ntype} matrix containing the occupation probability of each
#'  species on each patch type. Generally, each element is either a 0 or a 1, where 0 means the species cannot
#'  occupy that type and 1 means that it can.}
#'  \item{"type" - Patch types. A vector of integer with length equal to the number of patches. Each element 
#'  is an index to the patch types}
#' }
#' @return A vector of integers indexing the species occupying each patch after the model runs for one timestep
#' @export
stoch_patch <- function(t, y, params) {
  empt <- y == 0 #determine empty patches, 0 means empty, any other y means species with ID number 'y'
  #Add new mortality to empty patch vector; d=death rate
  dead <- c(sample(which(!empt), rbinom(1, sum(!empt), params$d)), which(empt)) 
  specnums <- table(y[!empt]) #count how many there are of each species
  specprobs2 <- rep(0, params$nspec) #initialize species probability vector
  nam <- as.numeric(names(specnums)) #grab the ID numbers of species
  #if (0 %in% nam){nums<-specnums[-which(nam==0)]}else{nums<-specnums} #Get rid of empty cells
  specprobs2[nam] <- specnums / sum(specnums) #Species prob is equal to proportion of cells occupied
  #add equal probability of being drawn from species pool. poolprob = prob of immigrant coming from pool
  specprobs <- (1 - params$poolprob) * specprobs2 + params$poolprob * (1 / params$nspec) 
  specmat <- matrix(specprobs, nrow = params$nspec, ncol = params$ntype) #put probs into species*resource matrix
  pmat <- specmat * params$occmat #multiply prob of immigrating by prob of being able to feed (0 or 1)
  pmatvec <- pmat[ , params$type[dead]] #pull out probs for all cells according to their resource type
  #Use probs to calculate which species should be placed into which empty cells if any
  y[dead] <- apply(pmatvec, 2, function(x) which(rmultinom(1, 1, x) == 1) * rbinom(1, 1, sum(x)))
  return(list(y))
}