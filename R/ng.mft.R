#' Multiple Fisher tests
#'
#' Perform Fisher's exact probability tests multiple times
#' @param cgt output from \code{ng.GOprep_TAIR()}, \code{[,1]:"locus", [,2]:"GOid"}
#' @param gn.test contig names for testing
#' @param alternative alternative hypothesis for the Fisher test
#' @return No. of genes, p.value, and calculation time
#' @examples
#' data(ulg)
#' head(ulg)
#' data(gl)
#' head(gl)
#' fisher.res <- ng.mft(ulg, gl)
#' head(fisher.res)
#' # The numbers of genes in the output table are given as follows.
#' #               gn.test
#' #             TRUE FALSE
#' # cgt   TRUE   xtt   xft   xnt
#' #      FALSE   xtf   xff   xnf
#' #              xtn   xfn   xnn
#' @author A.J. Nagano & Y. Sato
#' @export
#' @seealso [ng.GOprep_TAIR()]
ng.mft <- function(cgt,  gn.test, alternative="greater"){

  #cat(sprintf("%s\n", Sys.time()))

  gid.u <- unique(cgt[,"GOid"])

  ft.in <- matrix(0, nrow=length(gid.u), ncol=9)
  colnames(ft.in) <- c("xtt", "xft", "xtf", "xff", "xnt", "xnf", "xtn", "xfn", "xnn")
  rownames(ft.in) <- gid.u

  #               gn.test
  #             TRUE FALSE
  #Group  TRUE   xtt   xft   xnt
  #      FALSE   xtf   xff   xnf
  #              xtn   xfn   xnn

  ft.in[,"xnn"] <- length(unique(cgt[, "locus"]))

  gn.pp.gid <- table(cgt[, "GOid"])
  ft.in[names(gn.pp.gid), "xnt"] <- gn.pp.gid
  ft.in[,"xnf"] <- ft.in[,"xnn"] - ft.in[,"xnt"]

  ft.in[,"xtn"] <- length(intersect(gn.test, unique(cgt[, "locus"])))
  ft.in[,"xfn"] <- ft.in[,"xnn"] - ft.in[,"xtn"]

  gsea.test <- cgt[is.element(cgt[,"locus"], gn.test), ]
  gn.test.gid <- table(gsea.test[, "GOid"])
  ft.in[names(gn.test.gid), "xtt"] <- gn.test.gid

  ft.in[,"xtf"] <- ft.in[,"xtn"] - ft.in[,"xtt"]
  ft.in[,"xft"] <- ft.in[,"xnt"] - ft.in[,"xtt"]
  ft.in[,"xff"] <- ft.in[,"xnf"] - ft.in[,"xtf"]

  #cat(sprintf("%s\n", Sys.time()))

  #Fisher's exact test.  8? sec
  fr <- rep(1, nrow(ft.in))
  dt <- rep(1, nrow(ft.in))
  for(i in 1:nrow(ft.in)){
    start <- Sys.time()
    if(ft.in[i,"xtn"] > 1 && ft.in[i,"xnt"] > 1){
      contable <- matrix(ft.in[i, 1:4], ncol=2)
      tmp <- try(stats::fisher.test(contable, alternative=alternative), silent=TRUE)
      if(class(tmp)!="try-error") {
        fr[i] <- tmp$p.value
      } else {
        fr[i] <- NA
      }
    } else {
    }
    end <- Sys.time()
    dt[i] <- end - start
  }

  out <- cbind(fr, ft.in, dt)
  colnames(out) <- c("p.value", colnames(ft.in), "time")
  rownames(out) <- rownames(ft.in)

  #cat(sprintf("%s\n", Sys.time()))

  return(out)

}
