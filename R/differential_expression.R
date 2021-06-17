#' DE test using limma
#'
#' A wrapper for limma DE test.
#' @import edgeR
#' @import limma
#' @export
de.test.limma = function(counts, group, norm.factors = NULL) {
    y = edgeR::DGEList(counts = counts, group = group)
    if (is.null(norm.factors)) {
        y = edgeR::calcNormFactors(y)
    } else {
        y$samples$norm.factors = norm.factors
    }
    design = model.matrix(~ group)
    logCPM = edgeR::cpm(y, log = TRUE, prior.count = 1)
    fit = limma::lmFit(logCPM, design)
    fit = limma::eBayes(fit, trend = TRUE)
    return(fit)
}

#' DE test using edgeR
#'
#' A wrapper for edgeR to run DE test with custom sample normalization factors,
#' i.e. library sizes.
#' @import edgeR
#' @export
de.test.edgeR = function(counts, group, norm.factors = NULL) {
    y = edgeR::DGEList(counts = counts, group = group)
    if (is.null(norm.factors)) {
        y = edgeR::calcNormFactors(y)
    } else {
        y$samples$norm.factors = norm.factors
    }
    design = model.matrix(~ group)
    y = edgeR::estimateDisp(y, design)
    fit = edgeR::glmQLFit(y, design = design)
    edgeR::glmQLFTest(fit, coef=2:length(levels(group)))
}

#' DE test using DESeq2
#'
#' A wrapper for DESeq to run DE test with custom sample normalization factors,
#' i.e. library sizes.
#'
#' @import DESeq2
#' @export
de.test.DESeq = function(counts, group, norm.factors = NULL) {
    dds.x = DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                           colData = data.frame('group' = group),
                                           design = ~ group)

    normFactorMatrix = matrix(rep(norm.factors,ncol(counts)), nrow=nrow(counts),byrow = FALSE)
    if (is.null(norm.factors)) {
        dds.x = DESeq2::estimateSizeFactors(dds.x)
    } else {
        DESeq2::normalizationFactors(dds.x) = normFactorMatrix
    }
    dds.x = DESeq2::estimateDispersions(dds.x, quiet=TRUE) %>%
        DESeq2::nbinomWaldTest(quiet = TRUE)
    return(dds.x)
}

