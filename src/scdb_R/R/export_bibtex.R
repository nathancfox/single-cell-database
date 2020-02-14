#' Create BibTeX citations for the requested datasets.
#' 
#' Resolves the DOI for each requested dataset and exports
#' the resulting list of BibTeX citations to outfile. Does
#' not export duplicate citations.
#' 
#' @param uuids Character vector of UUIDs corresponding to
#'   the desired datasets.
#' @param outfile Character vector of length 1. The path to
#'   the file the citations should be written to.
#' @param append Logical vector of length 1. If TRUE, the
#'   outfile will be appended to. If FALSE, the outfile will
#'   be overwritten.
#' @return NULL
#' @export
export_bibtex <- function(uuids, outfile, append=FALSE) {
    error_uuids <- c()
    invalid_uuids <- c()
    bibtex_citations <- c()
    for (i in seq_along(uuids)) {
        results <- tryCatch(
            {
                uuid_to_row(uuids[i])
                list(skip=FALSE,
                     invalid=c())
            }, error=function(cond) {
                if (!(uuids[i] %in% invalid_uuids)) {
                    # invalid_uuids <- c(invalid_uuids, uuids[i])
                    invalid <- uuids[i]
                }
                # skip <- TRUE
                list(skip=TRUE,
                     invalid=invalid)
            }, warning=function(cond) {
                if (!(uuids[i] %in% invalid_uuids)) {
                    # invalid_uuids <- c(invalid_uuids, uuids[i])
                    invalid_uuids <- uuids[i]
                }
                # skip <- TRUE
                list(skip=TRUE,
                     invalid=invalid)
            }
        )
        invalid_uuids <- c(invalid_uuids, results$invalid)
        if (results$skip) {
            next
        }
        results <- tryCatch(
            {
                bt <- get_bibtex(uuids[i])
                list(bt=bt,
                     skip=FALSE,
                     error=c())
            }, error=function(cond) {
                if (!(uuids[i] %in% error_uuids)) {
                    # error_uuids <- c(error_uuids, uuids[i])
                    error <- uuids[i]
                }
                # skip <- TRUE
                list(bt=c(),
                     skip=TRUE,
                     error=error)
            }
        )
        error_uuids <- c(error_uuids, results$error)
        if (results$skip) {
            next
        }
        bt <- results$bt
        bibtex_citations <- c(bibtex_citations, bt)
    }
    if (length(invalid_uuids) > 0) {
        cat("ERROR: The following UUIDs are invalid:\n")
        for (i in seq_along(invalid_uuids)) {
            cat(paste("    ", invalid_uuids[i], ",\n",
                      sep=""))
        }
        cat("\n")
    }
    if (length(error_uuids) > 0) {
        cat("ERROR: The following UUIDs produced a non-200 HTTP status code:\n")
        for (i in seq_along(error_uuids)) {
            cat(paste("    ", error_uuids[i], ",\n",
                      sep=""))
        }
        cat("\n")
    }
    bibtex_citations <- unique(bibtex_citations)
    for (i in seq_along(bibtex_citations)) {
        if (i == 1) {
            cat(bibtex_citations[i], file=outfile, append=append)
        } else {
            cat(bibtex_citations[i], file=outfile, append=TRUE)
        }
        cat("\n\n", file=outfile, append=TRUE)
    }
    if (length(invalid_uuids) > 0 || length(error_uuids) > 0) {
        if (length(bibtex_citations) > 0) {
            cat("The remaining UUIDs were exported successfully")
        }
    }
}
