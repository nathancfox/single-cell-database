#' Get a loom filename.
#' 
#' \code{get_loom_filename()} returns the full filepath to the
#' loom file associated with the passed UUID.
#' 
#' This filepath is retrieved from the external metadata file.
#' If the file does not exist, an error will be thrown.
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return Character vector of length 1 containing the full path to
#'   the desired loom file.
get_loom_filename <- function(uuid) {
    df = read.table(get_PATH_TO_METADATA(),
                    header = TRUE,
                    sep = "\t",
                    quote = "",
                    row.names = NULL,
                    stringsAsFactors = FALSE)
    if (!(uuid %in% df$uuid)) {
        dir_entries <- list.dirs(get_PATH_TO_DATABASE())
        for (i in seq_along(dir_entries)) {
            if (dir_entries[i] == uuid
                    && 'expr_mat.loom' %in% list.files(file.path(get_PATH_TO_DATABASE(), dir_entries[i]))) {
                return(file.path(get_PATH_TO_DATABASE(), dir_entries[i], 'expr_mat.loom'))
            }
        }
        stop("uuid is not valid!")
    }
    filename <- df[df$uuid == uuid, "uuid"]
    filename = file.path(get_PATH_TO_DATABASE(), filename, "expr_mat.loom")
    if (!(file.exists(filename))) {
        stop("Retrieved filename does not exist!")
    }
    return(filename)
}
