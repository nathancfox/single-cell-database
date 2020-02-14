#' Get the expression matrix names.
#' 
#' Get the expression matrix names from a dataset. The first
#' entry will always be "matrix", referring to the expression
#' matrix stored under lfile[["matrix"]], not under
#' lfile[["layers"]].
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return A character vector containing the names of the matrices.
#' @export
get_expr_mat_names <- function(uuid) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    mat_names <- c("matrix", names(lfile[["layers"]]))
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(mat_names)
}
