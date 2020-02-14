#' Get the cell universal metadata.
#' 
#' Get the cell universal metadata from a given dataset
#' as a dataframe.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param keep_missing Logical vector of length 1. If TRUE, the returned
#'   dataframe will have all columns. If FALSE, columns where all values
#'   == -1 or "-1" (indicating a missing column) will be dropped from
#'   the returned dataframe.
#' @return A dataframe holding all the data available from the cell
#'   universal metadata.
#' @export
get_cell_univ <- function(uuid, keep_missing = TRUE) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    col_keys <- names(lfile[["col_attrs"]])
    # Initializes with the number of rows in the matrix,
    # even though the loom file stores the expression matrix
    # as genes x cells, and the SingleCellExperiment stores
    # cells as columns. This is because R stores data in 
    # Column-major order, also known as "Fortran-style",
    # but Python uses Row-major order, also known as "C-style".
    # This means that R will read all n-D datasets (where n > 1)
    # in a transposed orientation. Thus, although it should
    # initialize with the number of cells, it uses the number
    # of rows. This is to avoid having to run a transpose
    # operation just to get a dimension more readably.
    col_data <- data.frame(matrix(nrow = lfile[["matrix"]]$dims[1],
                                  ncol = 0))
    for (i in seq_along(col_keys)) {
        if (col_keys[i] == "CellID") {
            next
        }
        k <- paste("col_attrs/", col_keys[i], sep = "")
        if (length(lfile[[k]]$dims) == 1) {
            if (all(as.character(lfile[[k]][ ]) == "-1")) {
                if (keep_missing) {
                    col_data[col_keys[i]] <- lfile[[k]][ ]
                }
            } else {
                col_data[col_keys[i]] <- lfile[[k]][ ]
            }
        }
    }
    # HARDCODED CONSTANT FLAG
    column_order <- c("cluster",
                      "species",
                      "tissue",
                      "source_organism",
                      "sex",
                      "condition",
                      "batch",
                      "uuid")
    column_order <- column_order[column_order %in% colnames(col_data)]
    # If a dataframe is column_reordered this way with a single
    # column, it returns a vector instead of a single column
    # dataframe. Conveniently, if a dataframe only has a single
    # column, the columns are already sorted.
    if (length(column_order) > 1) {
        col_data <- col_data[, column_order]
    }
    rownames(col_data) <- lfile[["col_attrs/CellID"]][ ]
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(col_data)
}
