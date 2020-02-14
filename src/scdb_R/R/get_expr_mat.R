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
#' @importClassesFrom Matrix dgCMatrix
#' @export
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
