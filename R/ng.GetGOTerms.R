#' get GO terms
#'
#' Function to get GO term from GO ID
#' @param GOid input GO ID
#' @author A.J. Nagano & Y. Sato
#' @import GO.db AnnotationDbi
#' @export
ng.GetGOTerms <- function(GOid){

  out <- NULL
  for(i in GOid){
    tmp <- try(get(i, GO.db::GOTERM), silent=TRUE)
    if(class(tmp)=="try-error"){
      out <- c(out, "NA")
    } else {
      out <- c(out, AnnotationDbi::Term(tmp))
    }
  }
  return(out)
}
