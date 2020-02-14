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
#'   stored in the \code{S4Vectors::metadata()} slot, named "batch_key",
#'   "cell_author_annot", and "gene_author_annot".
#' @export
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
