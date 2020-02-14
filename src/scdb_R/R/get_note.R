#' Get the note for a dataset.
#' 
#' Retrieves the note for a dataset from the "notes.tsv"
#' file at the root of the database. These are typically
#' small notes used when an entry is still in progress or
#' being edited in the future.
#' 
#' @param uuid Character vector of length 1. The UUID of the
#'   desired dataset.
#' @return A character vector of length 1. If the UUID has
#'   a note, the note is returned. If not, an empty string
#'   is returned.
#' @export
get_note <- function(uuid) {
    notes = read.table(file.path(get_PATH_TO_DATABASE(), "notes.tsv"),
                       header=FALSE, sep="\t", stringsAsFactors=FALSE)
    colnames(notes) <- c("uuid", "note")
    rownames(notes) <- notes$uuid
    tryCatch(
        {
            get_loom_filename(uuid)
        },
        error=function(cond) {
            stop("Not a valid UUID!")
        }
    )
    if (uuid %in% notes$uuid) {
        return(notes[uuid, "note"])
    } else {
        return("")
    }
}
