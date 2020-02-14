#' Get the gene author-annotations.
#' 
#' Get the gene author-annotations from a given dataset
#' as a dataframe. Sets rownames preferentially as Accession,
#' Gene, or a number range in that order. Accession and Gene
#' will not be used if they have non-unique values.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return A dataframe holding all the data available from the gene
#'   author-annotations.
#' @export
get_gene_author_annot <- function(uuid) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    row_keys <- names(lfile[["gene_author_annot"]])
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
        k <- paste("gene_author_annot/", row_keys[i], sep = "")
        if (length(lfile[[k]]$dims) == 1) {
            row_data[row_keys[i]] <- lfile[[k]][ ]
        }
    }
    if (!hdf5r::h5attr(lfile[["row_attrs/Accession"]], "all_missing")) {
        rownames(row_data) <- lfile[["row_attrs/Accession"]][ ]
    } else if (!hdf5r::h5attr(lfile[["row_attrs/Gene"]], "all_missing")) {
        genes <- lfile[["row_attrs/Gene"]][ ]
        if (length(genes) == length(unique(genes))) {
            rownames(row_data) <- genes
        } else {
            # Do nothing and let the rownames be integers
        }
    } else {
        # Do nothing and let the rownames be integers
    }
    column_order <- hdf5r::h5attr(lfile[["gene_author_annot"]],
                                  "column_order")
    column_order <- strsplit(column_order, "|", fixed = TRUE)[[1]]
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
