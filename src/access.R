PATH_TO_METADATA <- "/data/single_cell_database/external_metadata.tsv"

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
    return(filename)
}

get_h5_conn <- function(uuid) {
    lfile <- hdf5r::H5File$new(get_loom_filename(uuid),
                               mode = "r")
    return(lfile)
}

get_loom_conn <- function(uuid) {
    # loomR fails to validate because it does not have
    # updated specs for loom 3.0, however if you
    # turn off the validation, it throws a warning.
    # So I'm suppressing the warning.
    suppressWarnings(
        lfile <- loomR::connect(get_loom_filename(uuid),
                                mode = "r",
                                skip.validate = TRUE)
    )
    return(lfile)
}

get_sce <- function(uuid,
                    assay_for_matrix = "counts",
                    counts_assay = NULL,
                    logcounts_assay = NULL) {
    lfile <- get_h5_conn(uuid)
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
    return(sce)
}