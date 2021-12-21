#' Preparing a list of GO terms with TAIR annotation files
#'
#' Making GO list from TAIR input file
#' @param fn path to the file name. The input file follows TAIR file format, named "ATH_GO_GOSLIM.txt".
#' @return A matrix listing AGI codes and GO terms
#' @author A.J. Nagano & Y. Sato
#' @import GO.db
#' @export
#' @seealso [ng.mft()]
ng.GOprep_TAIR = function(fn) {
  go <- utils::read.delim(fn, header=FALSE, as.is=TRUE)

  #make unique locus-GOid table
  #a <- aggregate(go[,6], by=list(go[,1]), FUN=unique)
  uagi <- unique(go[,1])
  ulg <- NULL
  progress <- 0
  start.time <- Sys.time()
  for(i in uagi){ #for 300 min
    tmp.goid <- go[go[,1]==i,6]
    tmp.goid <- unique(tmp.goid)

    go.out <- tmp.goid
    for(j in tmp.goid){
      r.try <- try(tmp <- get(j, GO.db::GOMFANCESTOR), silent=TRUE)
      if(class(r.try) != "try-error"){
        go.out <- c(go.out, tmp)
      }
      r.try <- try(tmp <- get(j, GO.db::GOCCANCESTOR), silent=TRUE)
      if(class(r.try) != "try-error"){
        go.out <- c(go.out, tmp)
      }
      r.try <- try(tmp <- get(j, GO.db::GOBPANCESTOR), silent=TRUE)
      if(class(r.try) != "try-error"){
        go.out <- c(go.out, tmp)
      }
    }
    go.out <- unique(go.out)
    tmp <- grep("GO:\\d{7}", go.out)
    go.out <- go.out[tmp]
    tmp <- cbind(rep(i, length(go.out)), go.out)
    ulg <- rbind(ulg, tmp)

    progress <- progress + 1
    cat(sprintf("%d/%d %s %s \n", progress, length(uagi), i, Sys.time()))
  }

  ulg <- unique(ulg)
  colnames(ulg) <- c("locus", "GOid")

  return(ulg)
}
