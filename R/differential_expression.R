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
