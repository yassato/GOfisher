#' Converting GOfisher outputs into REVIGO inputs
#'
#' This is an utility function to convert the \code{ng.mft()} output into the input of the "rrvigo" package.
#' @param r output table from \code{ng.mft()}
#' @param gn.test a set of genes to be examined
#' @param cgt output from \code{ng.GOprep_TAIR()}, \code{[,1]:"locus", [,2]:"GOid"}
#' @return A Data Frame object useful for the "rrvigo" package
#' @author Y. Sato
#' @examples
#' data(ulg)
#' data(gl)
#' fisher.res <- ng.mft(ulg, gl)
#' input <- ng.GOfisher2REVIGO(fisher.res[fisher.res[,"xtt"]>1,], gl, ulg)
#' @import GO.db
#' @seealso [ng.mft()]
#' @export
ng.GOfisher2REVIGO = function(r, gn.test, cgt) {
  pvalue <- (r[,"p.value"])
  qvalue <- stats::p.adjust(r[,"p.value"], method="BH")
  p.adj<-  stats::p.adjust(r[,"p.value"])

  tmp.id <- rownames(r)
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.gr <- paste0(r[,"xtt"],"/",r[,"xtn"])
  tmp.bg <- paste0(r[,"xnt"],"/",r[,"xnn"])
  tmp.geneID <- mapply(function(i) { return(paste(ng.SearchGOTerms(gn.test=gn.test,i,cgt=cgt),collapse="/")) }, tmp.id)
  names(tmp.geneID) = NULL

  out <- data.frame(tmp.id, tmp.description, tmp.gr, tmp.bg, pvalue, p.adj, qvalue, tmp.geneID ,r[,"xtt"])
  colnames(out) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")

  return(out)
}
