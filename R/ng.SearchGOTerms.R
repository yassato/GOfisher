#' search GO terms
#'
#' Function to search genes having a given GO ID
#' @param gn.test a set of genes to be examined
#' @param GOid query GO ID
#' @param cgt output from \code{ng.GOprep_TAIR()}, \code{[,1]:"locus", [,2]:"GOid"}
#' @author A.J. Nagano & Y. Sato
#' @examples
#' data(ulg)
#' data(gl)
#' res <- ng.SearchGOTerms(gl, "GO:0006952", ulg)
#'
#' data(des)
#' des[res,]
#' @export
ng.SearchGOTerms = function(gn.test, GOid, cgt) {
  withgo <- cgt[cgt[,"GOid"]==GOid, "locus"]
  gn.test = unique(gn.test)

  return(intersect(gn.test, withgo))
}
