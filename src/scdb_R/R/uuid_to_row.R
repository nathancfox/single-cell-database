#' Get an entry of the external metadata.
#' 
#' Get the entry of the external metadata corresponding
#' to the requested UUID.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param columns Character vector containing requested columns from
#'   the external metadata. If NULL, the entire entry will be returned.
#'   If a column is invalid, the function will still run, however the
#'   entry will contain "NOT A VALID COLUMN" in lieu of the actual value.
#'   Additionally, a warning will be thrown.
#' @return A dataframe holding one row (the entry) where the column
#'   names match the columns argument.
#' @export
uuid_to_row <- function(uuid, columns = NULL) {
    df <- get_extern_md()
    if (!(uuid %in% df$uuid)) {
        stop("UUID not found in metadata file!")
    }
    row <- df[df$uuid == uuid, ]
    if (is.null(columns)) {
        rownames(row) <- c(1)
        return(row)
    } else {
        values <- data.frame(matrix(nrow = 1, ncol = 0))
        warning_flag <- FALSE
        for (i in seq_along(columns)) {
            if (columns[i] %in% colnames(df)) {
                values[columns[i]] <- row[columns[i]]
                append(values, row[1, columns[i]])
            } else {
                if (!warning_flag) {
                    warning("At least one column is invalid!")
                }
                warning_flag <- TRUE
                values[columns[i]] <- "NOT A VALID COLUMN"
            }
        }
        return(values)
    }
}
