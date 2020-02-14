#' Get the external metadata.
#' 
#' Get the external metadata as a dataframe.
#' 
#' @return A dataframe containing the external metadata.
#' @export
get_extern_md <- function() {
    df <- read.table(get_PATH_TO_METADATA(),
                     header = TRUE,
                     sep = "\t",
                     quote = "",
                     stringsAsFactors = TRUE)
    df$title <- as.character(df$title)
    df$authors <- as.character(df$authors)
    df$abstract <- as.character(df$abstract)
    df$date_generated <- as.Date(df$date_generated)
    df$doi <- as.character(df$doi)
    df$accession <- as.character(df$accession)
    df$date_integrated <- as.Date(df$date_integrated)
    df$uuid <- as.character(df$uuid)
    return(df)
}
