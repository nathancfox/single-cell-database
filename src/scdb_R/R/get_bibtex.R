#' Get the BibTeX citation for an entry.
#' 
#' Resolves the DOI for a dataset to generate a BibTeX citation
#' as a string.
#' 
#' @param uuid Character vector of length 1. The UUID of the
#'   desired dataset.
#' @return A character vector of length 1. It contains the BibTeX
#'   entry as a printable string.
#' @export
get_bibtex <- function(uuid) {
    doi <- uuid_to_row(uuid)$doi
    response <- httr::GET(paste("https://doi.org/", doi, sep=""),
                          httr::add_headers('Accept' = 'application/x-bibtex'))
    if (response$status_code != 200) {
        stop(paste("Status Code is ", response$status_code,
                   ", not 200!", sep=""))
    }
    return(httr::content(response, "text", encoding="UTF-8"))
}
