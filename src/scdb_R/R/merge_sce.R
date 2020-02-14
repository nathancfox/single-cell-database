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
#' @export
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
