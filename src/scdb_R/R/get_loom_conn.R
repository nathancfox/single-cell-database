#' Get a loomR connection to a loom file.
#' 
#' \code{get_loom_conn()} returns an loomR
#' connection object to the loom file associated
#' with the passed UUID.
#' 
#' The only mode is read-only.
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param warning Logical vector of length 1. If \code{TRUE}, a long
#'   warning will be displayed, reminding the user that R accesses
#'   the loom file in column-major order, but that the file is stored
#'   in row-major order. Thus all accessed matrices will be transposed.
#' @return loomR connection to the desired loom file.
#' @export
get_loom_conn <- function(uuid, warning = TRUE) {
    # loomR fails to validate because it does not have
    # updated specs for loom 3.0, however if you
    # turn off the validation, it throws a warning.
    # So I'm suppressing the warning.
    suppressWarnings(
        lfile <- loomR::connect(get_loom_filename(uuid),
                                mode = "r",
                                skip.validate = TRUE)
    )
    if (warning) {
        cat(ARRAY_WARNING, "\n")
    }
    return(lfile)
}
