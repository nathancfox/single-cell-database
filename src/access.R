# Couldn't figure out how to use sparse matrices without
# explicitly attaching this package
library(Matrix)

# All places where something is hardcoded in instead of taking
# from global_constants.py are flagged with a comment
# reading "HARDCODED CONSTANT FLAG"

# HARDCODED CONSTANT FLAG
get_PATH_TO_DATABASE <- function() {
    if (Sys.info()[["nodename"]] == "dactyl.cshl.edu") {
        return("/tyronedata/single_cell_database/database")
    } else {
        return("/data/single_cell_database/database")
    }
}
# HARDCODED CONSTANT FLAG
get_PATH_TO_METADATA <- function() {
    if (Sys.info()[["nodename"]] == "dactyl.cshl.edu") {
        return("/tyronedata/single_cell_database/database/external_metadata.tsv")
    } else {
        return("/data/single_cell_database/database/external_metadata.tsv")
    }
}
ARRAY_WARNING <- paste("\nWARNING: R and Python interpret binary matrices\n",
                       "in transposed ways. R uses column-major order and\n",
                       "Python uses row-major order. Because the database\n",
                       "entries are created in Python, this means that all\n",
                       "n-D HDF5 datasets (where n > 1) accessed manually\n",
                       "from an HDF5 or loom connection will be returned\n",
                       "transposed, i.e. cells x genes. However, cell-specific\n",
                       "internal universal metadata will still be in the\n",
                       "\"col_attrs\" HDF5 group and vice versa for genes.\n\n",
                       "To silence this warning, pass the function arg:\n",
                       "    \"warning = FALSE\"\n",
                       sep = "")

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

#' Get an hdf5r connection to a loom file.
#' 
#' \code{get_h5_conn()} returns an hdf5r H5File
#' connection object to the loom file associated
#' with the passed UUID.
#' 
#' The only mode is read-only.
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param warning Logical vector of length 1. If \code{TRUE}, a long
#'   warning will be displayed, reminding the user that R accesses
#'   the loom file in column-major order, but that the file is stored
#'   in row-major order. Thus all accessed matrices will be transposed.
#' @return H5File connection to the desired loom file.
get_h5_conn <- function(uuid, warning = TRUE) {
    lfile <- hdf5r::H5File$new(get_loom_filename(uuid),
                               mode = "r")
    if (warning) {
        cat(ARRAY_WARNING, "\n")
    }
    return(lfile)
}

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

#' Get a SingleCellExperiment object from a loom file.
#' 
#' \code{get_sce()} returns a SingleCellExperiment object
#' containing the entire loom file loaded completely into
#' memory.
#' 
#' The loom file associated with the UUID is the one
#' loaded into an SCE object. The other 3 function arguments
#' allows any configuration of SingleCellExperiment assays.
#' 
#' @section Assay assignment:
#' The loom file has a main matrix in the HDF5 dataset 'matrix'.
#' The rest of the matrices, if available, are under the
#' HDF5 group 'layers', as datasets of their own. SingleCellExperiment
#' objects typically have the raw counts stored under the assay "counts",
#' and the log-normalized counts stored under the assay "logcounts".
#' Other matrices can be named anything. The \code{assay_for_matrix}
#' argument should hold the name of the assay that the 'matrix' dataset
#' should be stored in. Then, the \code{counts_assay} and
#' \code{logcounts_assay} arguments can hold the names of HDF5 layers
#' that should be stored in "counts" or "logcounts" assays respectively.
#' For example: if the loom file has a raw reads matrix in 'matrix'
#' and a normalized matrix in 'norm' and a 3rd matrix in 'misc', the
#' appropriate function call would be:
#' \code{get_sce(UUID,
#'               assay_for_matrix = "counts",
#'               counts_assay = NULL,
#'               logcounts_assay = "norm")}
#' However, if the 3rd matrix was in 'matrix' and the raw reads and
#' normalized reads were in 'raw' and 'norm', the appropriate function
#' call would be:
#' \code{get_sce(UUID,
#'               assay_for_matrix = "misc",
#'               counts_assay = "raw",
#'               logcounts_assay = "norm")}
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param assay_for_matrix Character vector of length 1. The name of the
#'   assay that the 'matrix' dataset should be stored in.
#' @param counts_assay Character vector of length 1, or NULL. If not NULL,
#'   must be the name of an HDF5 dataset under 'layers' that will be
#'   stored under the "counts" assay of the returned SingleCellExperiment
#'   object. For obvious reasons, this argument cannot be not NULL if
#'   \code{assay_for_matrix} == "counts".
#' @param logcounts_assay Character vector of length 1, or NULL. If not NULL,
#'   must be the name of an HDF5 dataset under 'layers' that will be
#'   stored under the "logcounts" assay of the returned SingleCellExperiment
#'   object. For obvious reasons, this argument cannot be not NULL if
#'   \code{assay_for_matrix} == "logcounts".
#' @return SingleCellExperiment object with the assays as assigned and the
#'   'col_attrs' datasets consolidated into
#'   \code{SingleCellExperiment::colData()} and the 'row_attrs' dataset
#'   consolidated into \code{SingleCellExperiment::rowData()}.
#'   Additionally, the batch_key, and author annotated metadata are
#'   stored in the \code{SummarizedExperiment::metadata()} slot,
#'   named "batch_key", "cell_author_annot", and "gene_author_annot".
get_sce <- function(uuid,
                    assay_for_matrix = "counts",
                    counts_assay = NULL,
                    logcounts_assay = NULL) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    layers <- names(lfile[["layers"]])

    # Assays
    if ("matrix" %in% layers) {
        stop(paste("ERROR: One of the layers is named \"matrix\"!"))
    }
    if (assay_for_matrix %in% layers) {
        stop(paste("assay_for_matrix cannot be the name of one ",
                   "of the layers in the dataset:\n",
                   "[", paste(layers, collapse = ", "), "]",
                   sep = ""))
    }
    if (!is.null(counts_assay)) {
        if (assay_for_matrix == "counts") {
            stop(paste("assay_for_matrix cannot be \"counts\" ",
                       "if counts_assay is not NULL",
                       sep = ""))
        }
        if (!(counts_assay %in% layers)) {
            stop(paste("counts_assay must be one of the layers ",
                    "in the dataset",
                    sep = ""))
        }
        if (counts_assay != "counts" && "counts" %in% layers) {
            stop(paste("counts_assay must be NULL if there is ",
                       "a layer named \"counts\" and it is not ",
                       "being assigned to the counts_assay"))
        }
    }
    if (!is.null(logcounts_assay)) {
        if (assay_for_matrix == "logcounts") {
            stop(paste("assay_for_matrix cannot be \"logcounts\" ",
                       "if logcounts_assay is not NULL",
                       sep = ""))
        }
        if (!(logcounts_assay %in% layers)) {
            stop(paste("logcounts_assay must be one of the layers ",
                    "in the dataset",
                    sep = ""))
        }
        if (logcounts_assay != "logcounts" && "logcounts" %in% layers) {
            stop(paste("logcounts_assay must be NULL if there is ",
                       "a layer named \"logcounts\" and it is not ",
                       "being assigned to the logcounts_assay"))
        }
    }
    if (!is.null(counts_assay) && !is.null(logcounts_assay)) {
        if (counts_assay == logcounts_assay) {
            stop("counts_assay and logcounts_assay cannot be equal")
        }
    }
    assays = list()
    assays[[assay_for_matrix]] <- t(lfile[["matrix"]][ , ])
    for (i in seq_along(layers)) {
        k <- paste("layers/", layers[i], sep = "")
        if (!is.null(counts_assay)) {
            if (layers[i] == counts_assay) {
                assays[["counts"]] <- t(lfile[[k]][ , ])
                next
            }
        } else if (!is.null(logcounts_assay)) {
            if (layers[i] == logcounts_assay) {
                assays[["logcounts"]] <- t(lfile[[k]][ , ])
                next
            }
        } else{
            assays[[layers[i]]] <- t(lfile[[k]][ , ])
        }
    }

    # Row Data
    row_data <- get_gene_univ(uuid)

    # Col Data
    col_data <- get_cell_univ(uuid)

    # Miscellaneous Metadata
    metadata <- list()
    metadata[["batch_key"]] <- get_batch_key(uuid)
    metadata[["cell_author_annot"]] <- get_cell_author_annot(uuid)
    metadata[["gene_author_annot"]] <- get_gene_author_annot(uuid)

    # Final Construction
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = assays,
        rowData = row_data,
        colData = col_data,
        metadata = metadata
    )
    lfile <- get_h5_conn(uuid, warning = FALSE)
    if (!hdf5r::h5attr(lfile[["row_attrs/Accession"]], "all_missing")) {
        rownames(sce) <- lfile[["row_attrs/Accession"]][ ]
    } else if (!hdf5r::h5attr(lfile[["row_attrs/Gene"]], "all_missing")) {
        genes <- lfile[["row_attrs/Gene"]][ ]
        if (length(genes) == length(unique(genes))) {
            rownames(sce) <- genes
        } else {
            # Do nothing and let the rownames be integers
        }
    } else {
        # Do nothing and let the rownames be integers
    }
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(sce)
}

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

#' Get the cell author-annotations.
#' 
#' Get the cell author-annotations from a given dataset
#' as a dataframe.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return A dataframe holding all the data available from the cell
#'   author-annotations.
get_cell_author_annot <- function(uuid) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    col_keys <- names(lfile[["cell_author_annot"]])
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
        k <- paste("cell_author_annot/", col_keys[i], sep = "")
        if (length(lfile[[k]]$dims) == 1) {
            col_data[col_keys[i]] <- lfile[[k]][ ]
        }
    }
    rownames(col_data) <- lfile[["col_attrs/CellID"]][ ]
    column_order <- hdf5r::h5attr(lfile[["cell_author_annot"]],
                                  "column_order")
    column_order <- strsplit(column_order, "|", fixed = TRUE)[[1]]
    # If a dataframe is column_reordered this way with a single
    # column, it returns a vector instead of a single column
    # dataframe. Conveniently, if a dataframe only has a single
    # column, the columns are already sorted.
    if (length(column_order) > 1) {
        col_data <- col_data[, column_order]
    }
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(col_data)
}

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

#' Get the external metadata.
#' 
#' Get the external metadata as a dataframe.
#' 
#' @return A dataframe containing the external metadata.
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

#' Get the cell IDs for a given dataset.
#' 
#' For a given dataset, retrieve the cell IDs.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return A character vector containing the cell IDs.
get_cell_ids <- function(uuid) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    cell_ids <- lfile[["col_attrs/CellID"]][ ]
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(cell_ids)
}

#' Get the gene IDs for a given dataset.
#' 
#' For a given dataset, retrieve the gene IDs.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param accession Logical vector of length 1. If TRUE, the "Accession"
#'   variable will be returned. If FALSE, the "Gene" variable will be
#'   returned. If the requested variable is not available, but the other
#'   is, the other will be returned with a warning.
#' @return A character vector containing the gene IDs.
get_gene_ids <- function(uuid, accession = TRUE) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    if (accession) {
        if (!hdf5r::h5attr(lfile[["row_attrs/Accession"]], "all_missing")) {
            gene_ids <- lfile[["row_attrs/Accession"]][ ]
        } else if (!hdf5r::h5attr(lfile[["row_attrs/Gene"]], "all_missing")) {
            warning("\"Accession\" not available. Returning \"Gene\" instead.")
            gene_ids <- lfile[["row_attrs/Gene"]][ ]
        } else {
            stop(paste("Dataset ", uuid, " does not have a \"Gene\"",
                       " or an \"Accession\" row attribute!", sep = ""))
        }
    } else {
        if (!hdf5r::h5attr(lfile[["row_attrs/Gene"]], "all_missing")) {
            gene_ids <- lfile[["row_attrs/Gene"]][ ]
        } else if (!hdf5r::h5attr(lfile[["row_attrs/Accession"]], "all_missing")) {
            warning("\"Gene\" not available. Returning \"Accession\" instead.")
            gene_ids <- lfile[["row_attrs/Accession"]][ ]
        } else {
            stop(paste("Dataset ", uuid, " does not have a \"Gene\"",
                       " or an \"Accession\" row attribute!", sep = ""))
        }
    }
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(gene_ids)
}

#' Get the description attribute for a given column.
#' 
#' For a given internal metadata column, retrieve the HDF5
#' "description" attribute, if available, from the requested
#' dataset.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param column Character vector of length 1. The name of the column, for
#'   which the "description" attribute should be retrieved.
#' @param var Character vector of length 1. Must be "cell" or "gene".
#'   Indicates which variable's metadata the column is in. 
#' @param metadata Character vector of length 1. Must be "universal" or
#'   "author_annot". Indicates which metadata the column is in.
#' @return A character vector of length 1, holding the HDF5 "description"
#'   attribute of the indicated column. If there is no description
#'   available, a message is printed, and the function returns NULL.
get_column_description <- function(uuid, column,
                                   var = "cell", metadata = "universal") {
    if (metadata == "universal") {
        if (var == "cell") {
            lfile <- get_h5_conn(uuid, warning = FALSE)
            if (column %in% names(lfile[["col_attrs"]])) {
                key <- paste("col_attrs/", column, sep = "")
                if ("description" %in% hdf5r::h5attr_names(lfile[[key]])) {
                    desc <- hdf5r::h5attr(lfile[[key]], "description")
                    if (is.null(desc) || desc == "") {
                        cat("No description available!\n")
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return()
                    } else {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return(desc)
                    }
                } else {
                    cat("No description available!\n")
                    # In hdf5r, closing a file closes all other file connection
                    # objects for that file, regardless of scope. This ensures
                    # that it's not left hanging, but allows for the possibility
                    # that it was closed elsewhere.
                    if (lfile$is_valid) {
                        lfile$close_all()
                    }
                    return()
                }
            } else {
                # In hdf5r, closing a file closes all other file connection
                # objects for that file, regardless of scope. This ensures
                # that it's not left hanging, but allows for the possibility
                # that it was closed elsewhere.
                if (lfile$is_valid) {
                    lfile$close_all()
                }
                stop(paste(column, " is not a valid column!", sep = ""))
            }
        } else if ( var == "gene") {
            lfile <- get_h5_conn(uuid, warning = FALSE)
            if (column %in% names(lfile[["row_attrs"]])) {
                key <- paste("row_attrs/", column, sep = "")
                if ("description" %in% hdf5r::h5attr_names(lfile[[key]])) {
                    desc <- hdf5r::h5attr(lfile[[key]], "description")
                    if (is.null(desc) || desc == "") {
                        cat("No description available!\n")
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return()
                    } else {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return(desc)
                    }
                } else {
                    cat("No description available!\n")
                    # In hdf5r, closing a file closes all other file connection
                    # objects for that file, regardless of scope. This ensures
                    # that it's not left hanging, but allows for the possibility
                    # that it was closed elsewhere.
                    if (lfile$is_valid) {
                        lfile$close_all()
                    }
                    return()
                }
            } else {
                # In hdf5r, closing a file closes all other file connection
                # objects for that file, regardless of scope. This ensures
                # that it's not left hanging, but allows for the possibility
                # that it was closed elsewhere.
                if (lfile$is_valid) {
                    lfile$close_all()
                }
                stop(paste(column, " is not a valid column!", sep = ""))
            }
        } else {
            stop("ERROR: var must be \"cell\" or \"gene\"!")
        }
    } else if(metadata == "author_annot") {
        if (var == "cell") {
            lfile <- get_h5_conn(uuid, warning = FALSE)
            if (column %in% names(lfile[["cell_author_annot"]])) {
                key <- paste("cell_author_annot/", column, sep = "")
                if ("description" %in% hdf5r::h5attr_names(lfile[[key]])) {
                    desc <- hdf5r::h5attr(lfile[[key]], "description")
                    if (is.null(desc) || desc == "") {
                        cat("No description available!\n")
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return()
                    } else {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return(desc)
                    }
                } else {
                    cat("No description available!\n")
                    # In hdf5r, closing a file closes all other file connection
                    # objects for that file, regardless of scope. This ensures
                    # that it's not left hanging, but allows for the possibility
                    # that it was closed elsewhere.
                    if (lfile$is_valid) {
                        lfile$close_all()
                    }
                    return()
                }
            } else {
                # In hdf5r, closing a file closes all other file connection
                # objects for that file, regardless of scope. This ensures
                # that it's not left hanging, but allows for the possibility
                # that it was closed elsewhere.
                if (lfile$is_valid) {
                    lfile$close_all()
                }
                stop(paste(column, " is not a valid column!", sep = ""))
            }
        } else if ( var == "gene") {
            lfile <- get_h5_conn(uuid, warning = FALSE)
            if (column %in% names(lfile[["gene_author_annot"]])) {
                key <- paste("gene_author_annot/", column, sep = "")
                if ("description" %in% hdf5r::h5attr_names(lfile[[key]])) {
                    desc <- hdf5r::h5attr(lfile[[key]], "description")
                    if (is.null(desc) || desc == "") {
                        cat("No description available!\n")
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return()
                    } else {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return(desc)
                    }
                } else {
                    cat("No description available!\n")
                    # In hdf5r, closing a file closes all other file connection
                    # objects for that file, regardless of scope. This ensures
                    # that it's not left hanging, but allows for the possibility
                    # that it was closed elsewhere.
                    if (lfile$is_valid) {
                        lfile$close_all()
                    }
                    return()
                }
            } else {
                # In hdf5r, closing a file closes all other file connection
                # objects for that file, regardless of scope. This ensures
                # that it's not left hanging, but allows for the possibility
                # that it was closed elsewhere.
                if (lfile$is_valid) {
                    lfile$close_all()
                }
                stop(paste(column, " is not a valid column!", sep = ""))
            }
        } else {
            stop("ERROR: var must be \"cell\" or \"gene\"!")
        }
    } else {
        stop("ERROR: metadata must be \"universal\" or \"author_annot\"!")
    }
}

#' Get an expression matrix from a dataset.
#' 
#' Get the named expression matrix from a dataset. Default
#' is to get the expression matrix stored at lfile[["matrix"]],
#' but any of the named layers can be gotten.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param matrix Character vector of length 1. Name of the desired expression
#'   matrix. If "matrix", the main expression matrix under lfile[["matrix"]]
#'   will be returned. Otherwise, must be a valid layer.
#' @return A genes x cells unnamed sparse matrix holding the requested
#'   expression matrix. It is a Matrix::dgCMatrix class.
get_expr_mat <- function(uuid, matrix = "matrix") {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    if (matrix == "matrix") {
        return_mat <- as(t(lfile[["matrix"]][ , ]), "dgCMatrix")
    } else {
        if (!(matrix %in% names(lfile[["layers"]]))) {
            stop("matrix must be \"matrix\" or a valid layer name!")
        } else {
            key <- paste("layers/", matrix, sep = "")
            return_mat <- as(t(lfile[[key]][ , ]), "dgCMatrix")
        }
    }
    colnames(return_mat) <- lfile[["col_attrs/CellID"]][ ]
    row_keys <- names(lfile[["row_attrs"]])
    if ("Accession" %in% row_keys) {
        if (length(lfile[["row_attrs/Accession"]][ ])
                == length(unique(lfile[["row_attrs/Accession"]][ ]))) {
            rownames(return_mat) <- lfile[["row_attrs/Accession"]][ ]
        } else {
            warning('\"Accession\" has non-unique values!')
            # Do Nothing and let the rownames be numbers
        }
   } else if("Gene" %in% row_keys) {
        if (length(lfile[["row_attrs/Gene"]][ ])
                == length(unique(lfile[["row_attrs/Gene"]][ ]))) {
            rownames(return_mat) <- lfile[["row_attrs/Gene"]][ ]
        } else {
            # Do Nothing and let the rownames be numbers
        }
    } else {
        warning("\"Gene\" and \"Accession\" are both missing!")
        # Do Nothing and let the rownames be numbers
    }
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(return_mat)
}

#' Get the "batch_key" attribute.
#' 
#' Get the "batch_key" attribute from the cell-specific
#' internal universal metadata field "batch" for the
#' indicated dataset.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return A character vector of length 1 containing the batch_key
get_batch_key <- function(uuid) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    if (!("batch" %in% names(lfile[["col_attrs"]]))) {
        stop("ERROR: Cell universal metadata doesn\'t hvae a \"batch\" column!")
    }
    if (!("batch_key" %in% hdf5r::h5attr_names(lfile[["col_attrs/batch"]]))) {
        stop(paste("ERROR: Cell universal metadata \"batch\" column ",
                   "doesn\'t have a \"batch_key\" attributes!",
                   sep = ""))
    }
    batch_key <- hdf5r::h5attr(lfile[["col_attrs/batch"]], "batch_key")
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(batch_key)
}

#' Get the expression matrix names.
#' 
#' Get the expression matrix names from a dataset. The first
#' entry will always be "matrix", referring to the expression
#' matrix stored under lfile[["matrix"]], not under
#' lfile[["layers"]].
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return A character vector containing the names of the matrices.
get_expr_mat_names <- function(uuid) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    mat_names <- c("matrix", names(lfile[["layers"]]))
    # In hdf5r, closing a file closes all other file connection
    # objects for that file, regardless of scope. This ensures
    # that it's not left hanging, but allows for the possibility
    # that it was closed elsewhere.
    if (lfile$is_valid) {
        lfile$close_all()
    }
    return(mat_names)
}

#' Get the all_missing attribute for a given column.
#' 
#' For a given internal metadata column, retrieve the HDF5
#' "all_missing" attribute, if available, from the requested
#' dataset.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @param column Character vector of length 1. The name of the column, for
#'   which the "all_missing" attribute should be retrieved.
#' @param var Character vector of length 1. Must be "cell" or "gene".
#'   Indicates which variable's metadata the column is in. 
#' @param metadata Character vector of length 1. Must be "universal" or
#'   "author_annot". Indicates which metadata the column is in.
#' @return A character vector of length 1, holding the HDF5 "all_missing"
#'   attribute of the indicated column. If there is no all_missing
#'   attribute available, a message is printed, and the function returns NULL.
get_column_allmissing <- function(uuid, column,
                                  var = "cell", metadata = "universal") {
    if (metadata == "universal") {
        if (var == "cell") {
            lfile <- get_h5_conn(uuid, warning = FALSE)
            if (column %in% names(lfile[["col_attrs"]])) {
                key <- paste("col_attrs/", column, sep = "")
                if ("all_missing" %in% hdf5r::h5attr_names(lfile[[key]])) {
                    all_missing <- hdf5r::h5attr(lfile[[key]], "all_missing")
                    if (!is.logical(all_missing)) {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        stop("all_missing is not a boolean!")
                    } else {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return(all_missing)
                    }
                } else {
                    # In hdf5r, closing a file closes all other file connection
                    # objects for that file, regardless of scope. This ensures
                    # that it's not left hanging, but allows for the possibility
                    # that it was closed elsewhere.
                    if (lfile$is_valid) {
                        lfile$close_all()
                    }
                    stop("all_missing attribute does not exist!")
                }
            } else {
                # In hdf5r, closing a file closes all other file connection
                # objects for that file, regardless of scope. This ensures
                # that it's not left hanging, but allows for the possibility
                # that it was closed elsewhere.
                if (lfile$is_valid) {
                    lfile$close_all()
                }
                stop(paste(column, " is not a valid column!", sep = ""))
            }
        } else if ( var == "gene") {
            lfile <- get_h5_conn(uuid, warning = FALSE)
            if (column %in% names(lfile[["row_attrs"]])) {
                key <- paste("row_attrs/", column, sep = "")
                if ("all_missing" %in% hdf5r::h5attr_names(lfile[[key]])) {
                    all_missing <- hdf5r::h5attr(lfile[[key]], "all_missing")
                    if (!is.logical(all_missing)) {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        stop("all_missing is not a boolean!")
                    } else {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return(all_missing)
                    }
                } else {
                    # In hdf5r, closing a file closes all other file connection
                    # objects for that file, regardless of scope. This ensures
                    # that it's not left hanging, but allows for the possibility
                    # that it was closed elsewhere.
                    if (lfile$is_valid) {
                        lfile$close_all()
                    }
                    stop("all_missing attribute does not exist!")
                }
            } else {
                # In hdf5r, closing a file closes all other file connection
                # objects for that file, regardless of scope. This ensures
                # that it's not left hanging, but allows for the possibility
                # that it was closed elsewhere.
                if (lfile$is_valid) {
                    lfile$close_all()
                }
                stop(paste(column, " is not a valid column!", sep = ""))
            }
        } else {
            stop("ERROR: var must be \"cell\" or \"gene\"!")
        }
    } else if(metadata == "author_annot") {
        if (var == "cell") {
            lfile <- get_h5_conn(uuid, warning = FALSE)
            if (column %in% names(lfile[["cell_author_annot"]])) {
                key <- paste("cell_author_annot/", column, sep = "")
                if ("all_missing" %in% hdf5r::h5attr_names(lfile[[key]])) {
                    all_missing <- hdf5r::h5attr(lfile[[key]], "all_missing")
                    if (!is.logical(all_missing)) {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        stop("all_missing is not a boolean!")
                    } else {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return(all_missing)
                    }
                } else {
                    # In hdf5r, closing a file closes all other file connection
                    # objects for that file, regardless of scope. This ensures
                    # that it's not left hanging, but allows for the possibility
                    # that it was closed elsewhere.
                    if (lfile$is_valid) {
                        lfile$close_all()
                    }
                    stop("all_missing attribute does not exist!")
                }
            } else {
                # In hdf5r, closing a file closes all other file connection
                # objects for that file, regardless of scope. This ensures
                # that it's not left hanging, but allows for the possibility
                # that it was closed elsewhere.
                if (lfile$is_valid) {
                    lfile$close_all()
                }
                stop(paste(column, " is not a valid column!", sep = ""))
            }
        } else if ( var == "gene") {
            lfile <- get_h5_conn(uuid, warning = FALSE)
            if (column %in% names(lfile[["gene_author_annot"]])) {
                key <- paste("gene_author_annot/", column, sep = "")
                if ("all_missing" %in% hdf5r::h5attr_names(lfile[[key]])) {
                    all_missing <- hdf5r::h5attr(lfile[[key]], "all_missing")
                    if (!is.logical(all_missing)) {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        stop("all_missing is not a boolean!")
                    } else {
                        # In hdf5r, closing a file closes all other file connection
                        # objects for that file, regardless of scope. This ensures
                        # that it's not left hanging, but allows for the possibility
                        # that it was closed elsewhere.
                        if (lfile$is_valid) {
                            lfile$close_all()
                        }
                        return(all_missing)
                    }
                } else {
                    # In hdf5r, closing a file closes all other file connection
                    # objects for that file, regardless of scope. This ensures
                    # that it's not left hanging, but allows for the possibility
                    # that it was closed elsewhere.
                    if (lfile$is_valid) {
                        lfile$close_all()
                    }
                    stop("all_missing attribute does not exist!")
                }
            } else {
                # In hdf5r, closing a file closes all other file connection
                # objects for that file, regardless of scope. This ensures
                # that it's not left hanging, but allows for the possibility
                # that it was closed elsewhere.
                if (lfile$is_valid) {
                    lfile$close_all()
                }
                stop(paste(column, " is not a valid column!", sep = ""))
            }
        } else {
            stop("ERROR: var must be \"cell\" or \"gene\"!")
        }
    } else {
        stop("ERROR: metadata must be \"universal\" or \"author_annot\"!")
    }
}

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

#' Get the BibTeX citation for an entry.
#' 
#' Resolves the DOI for a dataset to generate a BibTeX citation
#' as a string.
#' 
#' @param uuid Character vector of length 1. The UUID of the
#'   desired dataset.
#' @return A character vector of length 1. It contains the BibTeX
#'   entry as a printable string.
get_bibtex <- function(uuid) {
    doi <- uuid_to_row(uuid)$doi
    response <- httr::GET(paste("https://doi.org/", doi, sep=""),
                          httr::add_headers('Accept' = 'application/x-bibtex'))
    if (response$status_code != 200) {
        stop(paste("Status Code is ", response$status_code,
                   ", not 200!", sep=""))
    }
    return(httr::content(response, "text", encoding="UTF-8"))
}

#' Create BibTeX citations for the requested datasets.
#' 
#' Resolves the DOI for each requested dataset and exports
#' the resulting list of BibTeX citations to outfile. Does
#' not export duplicate citations.
#' 
#' @param uuids Character vector of UUIDs corresponding to
#'   the desired datasets.
#' @param outfile Character vector of length 1. The path to
#'   the file the citations should be written to.
#' @param append Logical vector of length 1. If TRUE, the
#'   outfile will be appended to. If FALSE, the outfile will
#'   be overwritten.
#' @return NULL
export_bibtex <- function(uuids, outfile, append=FALSE) {
    error_uuids <- c()
    invalid_uuids <- c()
    bibtex_citations <- c()
    for (i in seq_along(uuids)) {
        results <- tryCatch(
            {
                uuid_to_row(uuids[i])
                list(skip=FALSE,
                     invalid=c())
            }, error=function(cond) {
                if (!(uuids[i] %in% invalid_uuids)) {
                    # invalid_uuids <- c(invalid_uuids, uuids[i])
                    invalid <- uuids[i]
                }
                # skip <- TRUE
                list(skip=TRUE,
                     invalid=invalid)
            }, warning=function(cond) {
                if (!(uuids[i] %in% invalid_uuids)) {
                    # invalid_uuids <- c(invalid_uuids, uuids[i])
                    invalid_uuids <- uuids[i]
                }
                # skip <- TRUE
                list(skip=TRUE,
                     invalid=invalid)
            }
        )
        invalid_uuids <- c(invalid_uuids, results$invalid)
        if (results$skip) {
            next
        }
        results <- tryCatch(
            {
                bt <- get_bibtex(uuids[i])
                list(bt=bt,
                     skip=FALSE,
                     error=c())
            }, error=function(cond) {
                if (!(uuids[i] %in% error_uuids)) {
                    # error_uuids <- c(error_uuids, uuids[i])
                    error <- uuids[i]
                }
                # skip <- TRUE
                list(bt=c(),
                     skip=TRUE,
                     error=error)
            }
        )
        error_uuids <- c(error_uuids, results$error)
        if (results$skip) {
            next
        }
        bt <- results$bt
        bibtex_citations <- c(bibtex_citations, bt)
    }
    if (length(invalid_uuids) > 0) {
        cat("ERROR: The following UUIDs are invalid:\n")
        for (i in seq_along(invalid_uuids)) {
            cat(paste("    ", invalid_uuids[i], ",\n",
                      sep=""))
        }
        cat("\n")
    }
    if (length(error_uuids) > 0) {
        cat("ERROR: The following UUIDs produced a non-200 HTTP status code:\n")
        for (i in seq_along(error_uuids)) {
            cat(paste("    ", error_uuids[i], ",\n",
                      sep=""))
        }
        cat("\n")
    }
    bibtex_citations <- unique(bibtex_citations)
    for (i in seq_along(bibtex_citations)) {
        if (i == 1) {
            cat(bibtex_citations[i], file=outfile, append=append)
        } else {
            cat(bibtex_citations[i], file=outfile, append=TRUE)
        }
        cat("\n\n", file=outfile, append=TRUE)
    }
    if (length(invalid_uuids) > 0 || length(error_uuids) > 0) {
        if (length(bibtex_citations) > 0) {
            cat("The remaining UUIDs were exported successfully")
        }
    }
}

#' Merge SingleCellExperiment objects returned from the database.
#' 
#' Attempts to merge 2 SingleCellExperiment objects returned from
#' get_sce(). Will reject the merge if any incompatibilities
#' exist.
#' 
#' @param sce1 SingleCellExperiment object. The first
#'   SingleCellExperiment object to be merged. Should be an
#'   object returned from \code{get_sce()}.
#' @param sce2 SingleCellExperiment object. The second
#'   SingleCellExperiment object to be merged. Should be an
#'   object returned from \code{get_sce()}.
#' @param keep_author_annot Logical vector of length 1. If TRUE,
#'   the internal author-annotated metadata from each SCE will
#'   be labeled with cell_id_prefix and saved in the metadata
#'   of the new SCE. Otherwise, they will be dropped.
#' @param min_common_genes Numeric vector of length 1. Must
#'   be an non-negative integer. The two SCE objects must
#'   have this many valid gene IDs in common. Otherwise, the
#'   merge will fail.
#' @param cell_id_prefix Can be a single integer from the
#'   following list \code{c(4, 8, 12, 16, 20, 24, 28, 32)} or
#'   a character vector of length 2 with 2 different prefixes
#'   for the cell labels from \code{sce1} and \code{sce2}
#'   respectively. If it is an integer, that number of digits
#'   from the UUIDs are used as a prefix (does not include the
#'   hyphens). If the number is too small to produce unique
#'   prefixes, it is increased by 4 until it works.
#' @return A new SingleCellExperiment object containing all
#'   cells from both objects (\code{sce1} first). It only
#'   retains genes that were guaranteed to be unique inside
#'   each SCE object and in common between the two objects.
#'   Cell IDs are prefixed with \code{cell_id_prefix}. The
#'   \code{colData()} from each is also merged. The names of
#'   the author-annotated internal metadata dataframe, if
#'   kept are prefixed with the \code{cell_id_prefix}. The
#'   function attempts to combine the batch_keys into a 
#'   vector, but may not do this correctly if sce1 or sce2
#'   are the product of a previous merge.
merge_sce <- function(sce1, sce2,
                      keep_author_annot = FALSE,
                      min_common_genes = 15000,
                      cell_id_prefix = 4) {
    # Parse and merge gene IDs
    row_data_1 <- SingleCellExperiment::rowData(sce1)
    row_data_2 <- SingleCellExperiment::rowData(sce2)
    acc_1 <- NULL
    gene_1 <- NULL
    acc_2 <- NULL
    gene_2 <- NULL
    # Check for gene IDs in sce1
    if (!(all(row_data_1$Accession == "-1"))) {
        acc_1 <- row_data_1$Accession
    }
    if (!(all(row_data_1$Gene == "-1"))) {
        gene_1 <- row_data_1$Gene
    }
    if (is.null(acc_1) && is.null(gene_1)) {
        stop("sce1 has neither Accession nor Gene!")
    }
    # Check for gene IDs in sce2
    if (!(all(row_data_2$Accession == "-1"))) {
        acc_2 <- row_data_2$Accession
    }
    if (!(all(row_data_2$Gene == "-1"))) {
        gene_2 <- row_data_2$Gene
    }
    if (is.null(acc_2) && is.null(gene_2)) {
        stop("sce2 has neither Accession nor Gene!")
    }
    gene_var <- NULL
    if (is.null(acc_1) || is.null(acc_2)) {
        if (is.null(gene_1) || is.null(gene_2)) {
            stop("sce1 and sce2 do not have a gene ID variable in common!")
        } else {
            # Use Gene if possible
            if (length(gene_1) != length(unique(gene_1))) {
                stop(paste0("One of the objects does not have Accession ",
                            "IDs and sce1 only has non-unique gene IDs! ",
                            "Can't automatically merge!"))
            }
            if (length(gene_2) != length(unique(gene_2))) {
                stop(paste0("One of the objects does not have Accession ",
                            "IDs and sce2 only has non-unique gene IDs! ",
                            "Can't automatically merge!"))
            }
            gene_ids <- intersect(gene_1, gene_2)
            gene_var <- "Gene"
            if (length(gene_ids) < min_common_genes) {
                stop(paste0("There are only ",
                            length(gene_ids),
                            " genes in common! Failed to meet the ",
                            "min_common_genes function argument restriction!"))
            }
        }
    } else {
        # Use Accession if possible, then try Gene if not
        if (length(acc_1) != length(unique(acc_1))) {
            stop(paste0("Assumption violated! Accession in sce1 has ",
                        "non-unique values!"))
        }
        if (length(acc_2) != length(unique(acc_2))) {
            stop(paste0("Assumption violated! Accession in sce2 has ",
                        "non-unique values!"))
        }
        gene_ids <- intersect(acc_1, acc_2)
        gene_var <- "Accession"
        if (length(gene_ids) < min_common_genes) {
            if (is.null(gene_1) || is.null(gene_2)) {
                stop(paste0("There are only ",
                            length(gene_ids),
                            " genes in common! Failed to meet the ",
                            "min_common_genes function argument restriction!"))
            } else {
                if (length(gene_1) != length(unique(gene_1))) {
                    stop(paste0("There are only ",
                                length(gene_ids),
                                " Accession IDs in common! Failed to ",
                                "meet the min_common_genes function argument ",
                                "restriction. Gene IDs are not an option ",
                                "because sce1 only has non-unique Gene IDs! "))
                }
                if (length(gene_2) != length(unique(gene_2))) {
                    stop(paste0("There are only ",
                                length(gene_ids),
                                " Accession IDs in common! Failed to ",
                                "meet the min_common_genes function argument ",
                                "restriction. Gene IDs are not an option ",
                                "because sce2 only has non-unique Gene IDs! "))
                }
                gene_ids <- intersect(gene_1, gene_2)
                gene_var <- "Gene"
                if (length(gene_ids) < min_common_genes) {
                    stop(paste0("There are only ",
                                length(gene_ids),
                                " genes in common! Failed to meet the ",
                                "min_common_genes function argument restriction!"))
                }
            }
        }
    }
    # Figure out if the other gene field is coming along
    if (gene_var == "Accession") {
        if (!(is.null(gene_1) && !(is.null(gene_2)))) {
            # Bring Gene if they agree
            other_gene_ids_1 <- gene_1[match(gene_ids, acc_1)]
            other_gene_ids_2 <- gene_2[match(gene_ids, acc_2)]
            if (all(other_gene_ids_1 == other_gene_ids_2)) {
                other_gene_ids <- other_gene_ids_1
            } else {
                other_gene_ids <- NULL
            }
        } else {
            other_gene_ids <- NULL
        }
    } else if (gene_var == "Gene") {
        if (!(is.null(acc_1) && !(is.null(acc_2)))) {
            # Bring Gene if they agree
            other_gene_ids_1 <- acc_1[match(gene_ids, gene_1)]
            other_gene_ids_2 <- acc_2[match(gene_ids, gene_2)]
            if (all(other_gene_ids_1 == other_gene_ids_2)) {
                other_gene_ids <- other_gene_ids_1
            } else {
                other_gene_ids <- NULL
            }
        } else {
            other_gene_ids <- NULL
        }
    } else {
        stop(paste0("Assertion error! gene_var should be \"Accession\" ",
                    "or \"Gene\"! It is ",
                    gene_var))
    }
    # Reorder rows and reconstruct gene universal internal metadata
    sce1 <- sce1[match(gene_ids, row_data_1[[gene_var]]), ]
    sce2 <- sce2[match(gene_ids, row_data_2[[gene_var]]), ]
    rownames(sce1) <- gene_ids 
    rownames(sce2) <- gene_ids 
    SummarizedExperiment::rowData(sce1) <- NULL
    SummarizedExperiment::rowData(sce2) <- NULL
    SummarizedExperiment::rowData(sce1)[gene_var] <- gene_ids
    SummarizedExperiment::rowData(sce2)[gene_var] <- gene_ids
    if (!(is.null(other_gene_ids))) {
        if (gene_var == "Accession") {
            SummarizedExperiment::rowData(sce1)["Gene"] <- other_gene_ids
            SummarizedExperiment::rowData(sce2)["Gene"] <- other_gene_ids
        } else if (gene_var == "Gene") {
            SummarizedExperiment::rowData(sce1)["Accession"] <- other_gene_ids
            SummarizedExperiment::rowData(sce2)["Accession"] <- other_gene_ids
        }
    } else {
        if (gene_var == "Accession") {
            SummarizedExperiment::rowData(sce1)["Gene"] <- "-1"
            SummarizedExperiment::rowData(sce2)["Gene"] <- "-1"
        } else if (gene_var == "Gene") {
            SummarizedExperiment::rowData(sce1)["Accession"] <- "-1"
            SummarizedExperiment::rowData(sce2)["Accession"] <- "-1"
        }
    }
    rm(row_data_1)
    rm(row_data_2)
    rm(acc_1)
    rm(acc_2)
    rm(gene_1)
    rm(gene_2)
    # Reconstruct CellIDs and reorder columns
    col_data_1 <- SingleCellExperiment::colData(sce1)
    col_data_2 <- SingleCellExperiment::colData(sce2)
    cell_ids_1 <- rownames(col_data_1)
    cell_ids_2 <- rownames(col_data_2)
    if (is.numeric(cell_id_prefix)) {
        loop = TRUE
        adjust_flag = FALSE
        while(loop) {
            if (length(cell_id_prefix) != 1) {
                stop("If cell_id_prefix is a number, it must be a scalar!")
            }
            allowed_numbers <- c(4, 8, 12, 16, 20, 24, 28, 32)
            if (!(cell_id_prefix %in% allowed_numbers)) {
                stop(paste0("If cell_id_prefix is a number, it must be a multiple ",
                            "of 4 in the closed interval [4, 32]!"))
            } 
            if (cell_id_prefix <= 8) {
                pad <- 0
            } else if (cell_id_prefix <= 12) {
                pad <- 1
            } else if (cell_id_prefix <= 16) {
                pad <- 2
            } else if (cell_id_prefix <= 20) {
                pad <- 3
            } else {
                pad <- 4
            }
            cell_prefix_1 <- substr(col_data_1[1, "uuid"], 1, (cell_id_prefix + pad))
            cell_prefix_2 <- substr(col_data_2[1, "uuid"], 1, (cell_id_prefix + pad))
            if (cell_prefix_1 == cell_prefix_2) {
                if (!(adjust_flag)) {
                    old_cell_prefix <- cell_id_prefix
                    adjust_flag <- TRUE
                }
                if (cell_id_prefix == 32) {
                    stop("uuid of first cell in sce1 and sce2 is the same!")
                } else {
                    cell_id_prefix <- cell_id_prefix + 4
                }
            } else {
                loop <- FALSE
            }
        }
        if (adjust_flag) {
            warning(paste0("cell_id_prefix was not long enough to distinguish ",
                           "the two UUIDs. cell_id_prefix was lengthened to ",
                           cell_id_prefix,
                           "!"))
        }
    } else if (is.character(cell_id_prefix)) {
        if ((length(cell_id_prefix) != 2)
            || (length(cell_id_prefix) != length(unique(cell_id_prefix)))) {
            stop(paste0("If cell_id_prefix is a character vector, it must ",
                        "contain exactly 2 unique labels!"))
        }
        cell_prefix_1 <- cell_id_prefix[1]
        cell_prefix_2 <- cell_id_prefix[2]
    }
    cell_prefix_1 <- paste0(cell_prefix_1, "__")
    cell_prefix_2 <- paste0(cell_prefix_2, "__")
    cell_ids_1 <- paste0(cell_prefix_1, cell_ids_1)
    cell_ids_2 <- paste0(cell_prefix_2, cell_ids_2)
    colnames(sce1) <- cell_ids_1
    colnames(sce2) <- cell_ids_2
    sce1 <- sce1[, cell_ids_1]
    sce2 <- sce2[, cell_ids_2]
    # Remove un-mergeable Assays
    assays_1 <- SummarizedExperiment::assayNames(sce1)
    assays_2 <- SummarizedExperiment::assayNames(sce2)
    assays_to_keep <- intersect(assays_1, assays_2)
    assays_to_keep <- assays_to_keep[!(assays_to_keep %in% c("OTHER"))]
    if (length(assays_to_keep) == 0) {
        stop(paste0("There are no named assays in common! ",
                    "Automatic merge not possible!"))
    }
    for (i in seq_along(assays_1)) {
        if (!(assays_1[i] %in% assays_to_keep)) {
            SummarizedExperiment::assay(sce1, match(assays_1[i], assays_1)) <- NULL
            if (assays_1[i] == "OTHER") {
                warning(paste0("Assays labeled \"OTHER\" are ",
                               "automatically dropped!"))
            } else {
                warning(paste0("Dropped unmergeable assay \"",
                            assays_1[i],
                            "\" from sce1!"))
            }
        }
    }
    for (i in seq_along(assays_2)) {
        if (!(assays_2[i] %in% assays_to_keep)) {
            SummarizedExperiment::assay(sce2, match(assays_2[i], assays_2)) <- NULL
            if (assays_2[i] == "OTHER") {
                warning(paste0("Assays labeled \"OTHER\" are ",
                               "automatically dropped!"))
            } else {
                warning(paste0("Dropped unmergeable assay \"",
                            assays_2[i],
                            "\" from sce2!"))
            }
        }
    }
    # Save author annotated internal metadata if applicable
    if (keep_author_annot) {
        author_annot <- list(S4Vectors::metadata(sce1)[["cell_author_annot"]],
                             S4Vectors::metadata(sce1)[["gene_author_annot"]],
                             S4Vectors::metadata(sce2)[["cell_author_annot"]],
                             S4Vectors::metadata(sce2)[["gene_author_annot"]])
        names(author_annot) <- c(paste0(cell_prefix_1, "cell_author_annot"),
                                 paste0(cell_prefix_1, "gene_author_annot"),
                                 paste0(cell_prefix_2, "cell_author_annot"),
                                 paste0(cell_prefix_2, "gene_author_annot"))
    }
    # Reconstruct batch_key
    batch_key_1 <- S4Vectors::metadata(sce1)[["batch_key"]]
    if (length(batch_key_1) != 1) {
        warning(paste0("sce1 has more than one batch_key and may be the ",
                       "result of a previous merge. batch_key merging may ",
                       "be incorrect or messy!"))
    }
    if (is.null(names(batch_key_1))) {
        names(batch_key_1) <- rep(substr(cell_prefix_1,
                                         1, nchar(cell_prefix_1) - 2),
                                  length(batch_key_1))
    } else {
        names(batch_key_1) <- paste(substr(cell_prefix_1, 
                                           1, nchar(cell_prefix_1) - 2),
                                    names(batch_key_1))
    }
    batch_key_2 <- S4Vectors::metadata(sce2)[["batch_key"]]
    if (length(batch_key_2) != 1) {
        warning(paste0("sce2 has more than one batch_key and may be the ",
                       "result of a previous merge. batch_key merging may ",
                       "be incorrect or messy!"))
    }
    if (is.null(names(batch_key_2))) {
        names(batch_key_2) <- rep(substr(cell_prefix_2,
                                         1, nchar(cell_prefix_2) - 2),
                                  length(batch_key_2))
    } else {
        names(batch_key_2) <- paste(substr(cell_prefix_2,
                                           1, nchar(cell_prefix_2) - 2),
                                    names(batch_key_2))
    }
    batch_key <- c(batch_key_1, batch_key_2)
    # Clear metadata before the merge
    S4Vectors::metadata(sce1) <- list()
    S4Vectors::metadata(sce2) <- list()
    
    merged_sce <- SingleCellExperiment::cbind(sce1, sce2)
    if (keep_author_annot) {
        S4Vectors::metadata(merged_sce) <- author_annot
    }
    S4Vectors::metadata(merged_sce)[["batch_key"]] <- batch_key
    return(merged_sce)
}