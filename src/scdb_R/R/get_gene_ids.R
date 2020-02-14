#' Get the gene IDs for a given dataset.
#' 
#' For a given dataset, retrieve the gene IDs.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param accession Logical vector of length 1. If TRUE, the "Accession"
#'   variable will be returned. If FALSE, the "Gene" variable will be
#'   returned. If the requested variable is not available, but the other
#'   is, the other will be returned with a warning.
#' @return A character vector containing the gene IDs.
#' @export
get_gene_ids <- function(uuid, accession = TRUE) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    if (accession) {
        if (!hdf5r::h5attr(lfile[["row_attrs/Accession"]], "all_missing")) {
            gene_ids <- lfile[["row_attrs/Accession"]][ ]
        } else if (!hdf5r::h5attr(lfile[["row_attrs/Gene"]], "all_missing")) {
            warning("\"Accession\" not available. Returning \"Gene\" instead.")
            gene_ids <- lfile[["row_attrs/Gene"]][ ]
        } else {
            stop(paste("Dataset ", uuid, " does not have a \"Gene\"",
                       " or an \"Accession\" row attribute!", sep = ""))
        }
    } else {
        if (!hdf5r::h5attr(lfile[["row_attrs/Gene"]], "all_missing")) {
            gene_ids <- lfile[["row_attrs/Gene"]][ ]
        } else if (!hdf5r::h5attr(lfile[["row_attrs/Accession"]], "all_missing")) {
            warning("\"Gene\" not available. Returning \"Accession\" instead.")
            gene_ids <- lfile[["row_attrs/Accession"]][ ]
        } else {
            stop(paste("Dataset ", uuid, " does not have a \"Gene\"",
                       " or an \"Accession\" row attribute!", sep = ""))
        }
    }
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(gene_ids)
}
