#' Get an hdf5r connection to a loom file.
#' 
#' \code{get_h5_conn()} returns an hdf5r H5File
#' connection object to the loom file associated
#' with the passed UUID.
#' 
#' The only mode is read-only.
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param warning Logical vector of length 1. If \code{TRUE}, a long
#'   warning will be displayed, reminding the user that R accesses
#'   the loom file in column-major order, but that the file is stored
#'   in row-major order. Thus all accessed matrices will be transposed.
#' @return H5File connection to the desired loom file.
#' @export
get_h5_conn <- function(uuid, warning = TRUE) {
    lfile <- hdf5r::H5File$new(get_loom_filename(uuid),
                               mode = "r")
    if (warning) {
        cat(ARRAY_WARNING, "\n")
    }
    return(lfile)
}
