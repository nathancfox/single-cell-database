PATH_TO_METADATA <- "/data/single_cell_database/external_metadata.tsv"
ARRAY_WARNING <- paste("\nWARNING!! R and Python interpret binary matrices\n",
                       "in transposed ways. R uses column-major order and\n",
                       "Python uses row-major order. Because the database\n",
                       "entries are created in Python, this means that all\n",
                       "n-D HDF5 datasets (where n > 1) accessed manually\n",
                       "from an HDF5 or loom connection will be returned\n",
                       "transposed, i.e. cells x genes. However, cell metadata\n",
                       "will still be in the \"colattrs\" HDF5 group and\n",
                       "vice versa for genes.\n\n",
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
    df = read.table(PATH_TO_METADATA,
                    header = TRUE,
                    sep = "\t",
                    quote = "",
                    row.names = NULL,
                    stringsAsFactors = FALSE)
    if (!(uuid %in% df$uuid)) {
        stop("uuid is not valid!")
    }
    filename <- df[df$uuid == uuid, "file_location"]
    filename = file.path(filename, "expr_mat.loom")
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
get_sce <- function(uuid,
                    assay_for_matrix = "counts",
                    counts_assay = NULL,
                    logcounts_assay = NULL) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    layers <- names(lfile[["layers"]])

    # Assays
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
        } else if (counts_assay == "matrix") {
            warning(paste("counts assay is being assigned the ",
                          "layer named matrix, not the root dataset ",
                          "named matrix",
                          sep = ""))
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
        } else if (logcounts_assay == "matrix") {
            warning(paste("logcounts assay is being assigned the ",
                          "layer named matrix, not the root dataset ",
                          "named matrix",
                          sep = ""))
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
        if (row_keys[i] == "author_annot") {
            next
        }
        k <- paste("row_attrs/", row_keys[i], sep = "")
        if (length(lfile[[k]]$dims) == 1) {
            row_data[row_keys[i]] <- lfile[[k]][ ]
        }
    }
    if ("Gene" %in% row_keys) {
        rownames(row_data) <- lfile[["row_attrs/Gene"]][ ]
    } else if ("Accession" %in% row_keys) {
        rownames(row_data) <- lfile[["row_attrs/Accession"]][ ]
    } else {
        # Do Nothing
    }
    if (ncol(row_data) == 0) {
        row_data <- NULL
    }

    # Col Data
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
        if (col_keys[i] == "author_annot") {
            next
        }
        k <- paste("col_attrs/", col_keys[i], sep = "")
        if (length(lfile[[k]]$dims) == 1) {
            col_data[col_keys[i]] <- lfile[[k]][ ]
        }
    }
    rownames(col_data) <- lfile[["col_attrs/CellID"]][ ]
    if (ncol(col_data) == 0) {
        col_data <- NULL
    }

    # Final Construction
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = assays,
        rowData = row_data,
        colData = col_data
    )
    lfile$close_all()
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
        if (col_keys[i] == "author_annot") {
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
    rownames(col_data) <- lfile[["col_attrs/CellID"]][ ]
    if (ncol(col_data) == 0) {
        col_data <- NULL
    }
    lfile$close_all()
    return(col_data)
}

#' Get the gene universal metadata.
#' 
#' Get the gene universal metadata from a given dataset
#' as a dataframe.
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
        if (col_keys[i] == "author_annot") {
            next
        }
        k <- paste("row_attrs/", row_keys[i], sep = "")
        if (length(lfile[[k]]$dims) == 1) {
            if (all(as.character(lfile[[k]][ ]) == "-1")) {
                if (keep_missing) {
                    row_data[row_keys[i]] <- lfile[[k]][ ]
                }
            } else {
                row_data[row_keys[i]] <- lfile[[k]][ ]
            }
        }
    }
    if ("Gene" %in% row_keys) {
        rownames(row_data) <- lfile[["row_attrs/Gene"]][ ]
    } else if ("Accession" %in% row_keys) {
        rownames(row_data) <- lfile[["row_attrs/Accession"]][ ]
    } else {
        # Do Nothing
    }
    if (ncol(row_data) == 0) {
        row_data <- NULL
    }
    lfile$close_all()
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
    col_keys <- names(lfile[["col_attrs/author_annot"]])
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
        k <- paste("col_attrs/author_annot/", col_keys[i], sep = "")
        if (length(lfile[[k]]$dims) == 1) {
            col_data[col_keys[i]] <- lfile[[k]][ ]
        }
    }
    rownames(col_data) <- lfile[["col_attrs/CellID"]][ ]
    if (ncol(col_data) == 0) {
        col_data <- NULL
    }
    lfile$close_all()
    return(col_data)
}

#' Get the gene author-annotations.
#' 
#' Get the gene author-annotations from a given dataset
#' as a dataframe.
#' 
#' @param uuid Character vector of length 1. UUID of the desired dataset.
#' @return A dataframe holding all the data available from the gene
#'   author-annotations.
get_gene_author_annot <- function(uuid) {
    lfile <- get_h5_conn(uuid, warning = FALSE)
    row_keys <- names(lfile[["row_attrs/author_annot"]])
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
        k <- paste("row_attrs/author_annot/", row_keys[i], sep = "")
        if (length(lfile[[k]]$dims) == 1) {
            row_data[row_keys[i]] <- lfile[[k]][ ]
        }
    }
    if ("Gene" %in% row_keys) {
        rownames(row_data) <- lfile[["row_attrs/Gene"]][ ]
    } else if ("Accession" %in% row_keys) {
        rownames(row_data) <- lfile[["row_attrs/Accession"]][ ]
    } else {
        # Do Nothing
    }
    if (ncol(row_data) == 0) {
        row_data <- NULL
    }
    lfile$close_all()
    return(row_data)
}

get_extern_md <- function() {
    df <- read.table(PATH_TO_METADATA,
                     header = TRUE,
                     sep = "\t",
                     quote = "",
                     stringsAsFactors = TRUE)
    df$condition <- as.logical(df$condition)
    df$date_generated <- as.Date(df$date_generated)
    df$umis <- as.logical(df$umis)
    df$spikeins <- as.logical(df$spikeins)
    df$doi <- as.character(df$doi)
    df$accession <- as.character(df$accession)
    df$date_integrated <- as.Date(df$date_integrated)
    df$uuid <- as.character(df$uuid)
    df$file_location <- as.character(df$file_location)
    df$internal <- as.logical(df$internal)
    return(df)
}