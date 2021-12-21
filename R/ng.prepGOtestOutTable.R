#' Preparing an output table for GO tests
#'
#' Result table for Fisher tests, with multiple testing corrected
#' @param r output table from \code{ng.mft()}
#' @param adj.method A method for the correction of multiple testing. Passed to \code{p.adjust()}
#' @param alpha cut-off value of the false discovery rate
#' @author A.J. Nagano & Y. Sato
#' @examples
#' data(ulg)
#' data(gl)
#' fisher.res <- ng.mft(ulg, gl)
#' GO.list <- ng.prepGOtestOutTable(fisher.res[fisher.res[,"xtt"]>1,], alpha=0.05)
#' GO.list[order(GO.list[,1]),]
#' @export
#' @seealso [ng.mft()] [p.adjust()]
ng.prepGOtestOutTable <- function(r, adj.method="BH", alpha=0.01){

  adp <- stats::p.adjust(r[,"p.value"], method=adj.method)

  tmp.id <- rownames(r)[adp < alpha]
  tmp.adp <- adp[adp < alpha]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- r[adp < alpha, c("xtt", "xtn", "xnt", "xnn")]

  out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  colnames(out) <- c("p-value", "ID", "Description",  "A & B", "A",	"B", "U")

  return(out)
}
