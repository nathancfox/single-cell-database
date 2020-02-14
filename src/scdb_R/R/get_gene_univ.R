#' Get the gene universal metadata.
#' 
#' Get the gene universal metadata from a given dataset
#' as a dataframe. SEts rownames preferentially as Accession,
#' Gene, or a number range in that order. Accession and Gene
#' will not be used if they have non-unique values.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param keep_missing Logical vector of length 1. If TRUE, the returned
#'   dataframe will have all columns. If FALSE, columns where all values
#'   == -1 or "-1" (indicating a missing column) will be dropped from
#'   the returned dataframe.
#' @return A dataframe holding all the data available from the gene
#'   universal metadata.
#' @export
get_gene_univ <- function(uuid, keep_missing = TRUE) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    row_keys <- names(lfile[["row_attrs"]])
    # Initializes with the number of columns in the matrix,
    # even though the loom file stores the expression matrix
    # as genes x cells, and the SingleCellExperiment stores
    # genes as rows. This is because R stores data in 
    # Column-major order, also known as "Fortran-style",
    # but Python uses Row-major order, also known as "C-style".
    # This means that R will read all n-D datasets (where n > 1)
    # in a transposed orientation. Thus, although it should
    # initialize with the number of genes, it uses the number
    # of columns. This is to avoid having to run a transpose
    # operation just to get a dimension more readably.
    row_data <- data.frame(matrix(nrow = lfile[["matrix"]]$dims[2],
                                  ncol = 0))
    for (i in seq_along(row_keys)) {
        k <- paste("row_attrs/", row_keys[i], sep = "")
        if (length(lfile[[k]]$dims) == 1) {
            if (hdf5r::h5attr(lfile[[k]], "all_missing")) {
                if (keep_missing) {
                    row_data[row_keys[i]] <- lfile[[k]][ ]
                }
            } else {
                row_data[, row_keys[i]] <- lfile[[k]][ ]
            }
        }
    }
    # HARDCODED CONSTANT FLAG
    column_order <- c("Accession",
                      "Gene")
    column_order <- column_order[column_order %in% colnames(row_data)]
    # If a dataframe is column_reordered this way with a single
    # column, it returns a vector instead of a single column
    # dataframe. Conveniently, if a dataframe only has a single
    # column, the columns are already sorted.
    if (length(column_order) > 1) {
        row_data <- row_data[, column_order]
    }
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(row_data)
}
