obtain_gene_info <- function(genelist, genepos) {

    if (!any(genelist %in% rownames(genepos))) {
        stop("None of the genes can be found the position information. Please check the data.")
    } else if (!all(genelist %in% rownames(genepos))) {
        genemiss <- genelist[!(genelist %in% rownames(genepos))]
        warning(sprintf("The following genes cannot be found the position information: %s.", 
                        paste(genemiss, collapse=", ")))
    }
    df <- genepos[rownames(genepos) %in% genelist, ]
    df <- df[order(match(rownames(df), genelist)), ]
    
    return (df)
    
}

match_gene_snp <- function(geneinfo, snppos, cisrange) {
    emptyrecord <- NULL
    for (i in 1:nrow(geneinfo)) {
        cispos <- which(snppos[, 1] == geneinfo[i, 1] &
                            snppos[, 2] > (geneinfo[i, 2] - cisrange) &
                            snppos[, 2] < (geneinfo[i, 2] + cisrange))
        if (length(cispos) == 0) {
            emptyrecord <- c(emptyrecord, i)
        }
    }
    if (length(emptyrecord) == nrow(geneinfo)) {
        stop("None of the genes can be matched with SNPs. Please check the data.")
    } else if (length(emptyrecord) > 0) {
        warning(sprintf("The following genes cannot be matched with SNPs: %s.", 
                        paste(rownames(geneinfo)[emptyrecord], collapse=", ")))
        geneinfo <- geneinfo[-emptyrecord, ]
    }
    return (geneinfo)
}



