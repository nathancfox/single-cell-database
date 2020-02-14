#' Gets the path to the external metadata.
#'
#' Gets the full filepath to the external metadata file.
#' Takes into account the different naming on dactyl.
#'
#' @return A character vector of length 1 holding the full
#'   path to the external metadata file.
get_PATH_TO_METADATA <- function() {
# HARDCODED CONSTANT FLAG
    if (Sys.info()[["nodename"]] == "dactyl.cshl.edu") {
        return("/tyronedata/single_cell_database/database/external_metadata.tsv")
    } else {
        return("/data/single_cell_database/database/external_metadata.tsv")
    }
}
