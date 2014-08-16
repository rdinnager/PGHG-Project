#' Function to take a vector of species IDs and species names, and calculate the number of each species 
#' occuring in the vector. For use with \code{\link{stoch_patch}} model output.
#' @param spec_vec An integer vector of species IDs
#' @param spec_names A character vector of species names to which the \code{spec_vec} acts as an index
#' @return A named vector, tabulating the number of times each species occurs in the ID vector.
#' @export 
c_tabulate <- function(spec_vec, spec_names) {
  lev <- as.character(seq_along(spec_names))
  spec_tab <- table(factor(spec_vec, levels = lev))
  names(spec_tab) <- spec_names
  return(spec_tab)
}
