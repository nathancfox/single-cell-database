#' Get the "batch_key" attribute.
#' 
#' Get the "batch_key" attribute from the cell-specific
#' internal universal metadata field "batch" for the
#' indicated dataset.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return A character vector of length 1 containing the batch_key
#' @export
get_batch_key <- function(uuid) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    if (!("batch" %in% names(lfile[["col_attrs"]]))) {
        stop("ERROR: Cell universal metadata doesn\'t hvae a \"batch\" column!")
    }
    if (!("batch_key" %in% hdf5r::h5attr_names(lfile[["col_attrs/batch"]]))) {
        stop(paste("ERROR: Cell universal metadata \"batch\" column ",
                   "doesn\'t have a \"batch_key\" attributes!",
                   sep = ""))
    }
    batch_key <- hdf5r::h5attr(lfile[["col_attrs/batch"]], "batch_key")
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(batch_key)
}
