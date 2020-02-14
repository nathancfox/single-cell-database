#' Get the path to the root of the database
#' 
#' Gets the full file path to the root of the database.
#' Takes into account the differential naming on dactyl.
#'
#' @return A character vector of length 1 holding the full
#'   path to the root directory of the database.
get_PATH_TO_DATABASE <- function() {
# HARDCODED CONSTANT FLAG
    if (Sys.info()[["nodename"]] == "dactyl.cshl.edu") {
        return("/tyronedata/single_cell_database/database")
    } else {
        return("/data/single_cell_database/database")
    }
}
