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
#' @export
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
