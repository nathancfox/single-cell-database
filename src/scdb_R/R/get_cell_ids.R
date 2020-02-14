#' Get the cell IDs for a given dataset.
#' 
#' For a given dataset, retrieve the cell IDs.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return A character vector containing the cell IDs.
#' @export
get_cell_ids <- function(uuid) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    cell_ids <- lfile[["col_attrs/CellID"]][ ]
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(cell_ids)
}
