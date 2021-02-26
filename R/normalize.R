
#' normalize.by.qsmooth
#'
#' @import qsmooth
#' @export
normalize.by.qsmooth <- function(X, group_factor, ...) {
    qs <- qsmooth::qsmooth(t(X), group_factor = group_factor, ...)
    qs@qsmoothData %>%
        t() %>%
        return()
}

#' normalize.by.tmm
#'
#' Call calcNormFactors by edgeR
#' @import edgeR
#' @param X read count matrix with samples in rows and genes in columns
#' @export
normalize.by.tmm <- function(X,...) {
    normfactors.tmm = apply(X, 1, sum) / edgeR::calcNormFactors(t(X), method = 'TMM',...) # effective library sizes
    return(normalize.by.scaleFactors(X, normfactors.tmm))
    # X.normed = t(sapply(1:length(normfactors.tmm), FUN = function(j) {return(X[j,]/normfactors.tmm[j])}))
    # rownames(X.normed) = rownames(X)
    # return(X.normed)
}

#' normalize.by.deseq
#'
#' Call estimateSizeFactorsForMatrix by DESeq2
#' @import DESeq2
#' @param X read count matrix in the form genes x samples
#' @export
normalize.by.deseq <- function(X, ...) {
    normFactors = DESeq2::estimateSizeFactorsForMatrix(X)
    return(normalize.by.scaleFactors(X, normFactors))
    # X.normed = t(sapply(1:length(normFactors), FUN = function(j) {return(X[j,]/normFactors[j])}))
    # rownames(X.normed) = rownames(X)
    # return(X.normed)
}

#' normalize.by.poissonseq
#'
#' @import PoissonSeq
#' @export
normalize.by.poissonseq <- function(X, ...) {
    PoissonSeq::PS.Est.Depth(t(X), ...) %>%
        normalize.by.scaleFactors(X, .) %>%
        return()
}

#' normalize.by.refs
#'
#' normalize a read count matrix given the set of reference genes identified by id
#' @param X a read-count matrix of the form genes x samples
#' @param ref.idx an integer vector specifying the column indices of X to be used as reference
#' @export
normalize.by.refs <- function(X, ref.idx, scale = TRUE) {
    if (length(ref.idx) == 1) { # need to be treated specially since dim(X) reduces to NULL and cause error in apply
        Xref = matrix(X[ref.idx,], nrow=1)
    } else {
        Xref = X[ref.idx,]
    }
    normFactors = colSums(Xref) # apply(Xref,MARGIN = 2,FUN = sum)
    if (scale) {
        normFactors = normFactors / geom.mean(normFactors)
    }
    # sanity check
    idx.zero = which(normFactors == 0)
    if (length(idx.zero) > 0) {
        message(paste0("All reference transcripts are zero in the sample ", idx.zero, ". Please remove.\n"))
        return(NULL)
    }

    # X.norm = sapply(1:length(normFactors), FUN = function(i) {
    #     return(X[,i] / normFactors[i])
    # })
    X.norm = sweep(X, MARGIN = 2, STATS = normFactors, FUN = '/')
    colnames(X.norm) = colnames(X)
    rownames(X.norm) = rownames(X)
    return(X.norm)
}

#' normalize.by.scaleFactors
#'
#' @param X count matrix in the form genes x samples
#' @scaleFactors a vector of scaling factor, must be the same length as the number of samples
#' @export
normalize.by.scaleFactors <- function(X, scaleFactors) {
    if (length(scaleFactors) != ncol(X)) {
        stop("scaleFactors should have the same length as number of samples in the input matrix.")
    }
    X.normed = sapply(1:length(scaleFactors), FUN = function(j) {
        return(X[,j]/scaleFactors[j])
        })
    `colnames<-`(X.normed, colnames(X))
    return(X.normed)
}

uq.v1 <- function(X, group, ...) {
    effLibSizes = colSums(X) * edgeR::calcNormFactors(X, method = 'upperquartile', group = group, ...) # effective library sizes
    sweep(X, 2, effLibSizes, "/") %>%
        return()
}
uq.v2 <- function(X, group, ...) {
    effLibSizes = colSums(X) * edgeR::calcNormFactors(X, method = 'upperquartile', group = group, ...) # effective library sizes
    sweep(X, 2, mean(effLibSizes) / effLibSizes, "*") %>%
        return()
}
ruv_r.1 <- function(X, group) {
    if (!is.factor(group)) group <- factor(group)
    design <- model.matrix(~group, data=as.data.frame(X))

    y <- edgeR::DGEList(counts=X, group=group)
    y <- edgeR::calcNormFactors(y, method="upperquartile")
    y <- edgeR::estimateGLMCommonDisp(y, design)
    y <- edgeR::estimateGLMTagwiseDisp(y, design)

    fit <- edgeR::glmFit(y, design)
    res <- residuals(fit, type="deviance")
    if (is.null(res)) {
        message(str(y))
    }
    X.normed <- RUVSeq::RUVr(X, 1:nrow(X), k=1, res, round = FALSE)$normalizedCounts
    return(X.normed)
}

tmm.v1 <- function(X, group, ...) {
    effLibSizes = colSums(X) * edgeR::calcNormFactors(X, method = 'TMM', group = group, ...) # effective library sizes
    sweep(X, 2, effLibSizes, "/") %>%
        return()
}

tmm.v2 <- function(X, group, ...) {
    effLibSizes = colSums(X) * edgeR::calcNormFactors(X, method = 'TMM', group = group, ...) # effective library sizes
    sweep(X, 2, mean(effLibSizes) / effLibSizes, "*") %>%
        return()
}

pseq.v1 <- function(X, group = NA, ...) {
    PoissonSeq::PS.Est.Depth(X, ...) %>%
        sweep(X, 2, ., "/") %>%
        return()
}

deseq.v1 <- function(X, group) {
    DESeq2::estimateSizeFactorsForMatrix(X) %>%
        sweep(X, 2, ., "/")  %>%
        return()
}

deges.3 <- function(X, group, norm.method = 'tmm', test.method = 'edger', iteration = 1) {
    tcc <- TCC::TCC(X, group) %>%
        TCC::calcNormFactors(norm.method = 'tmm', test.method = 'edger', iteration = iteration) %>%
        TCC::getNormalizedData() %>%
        return()
}

tc <- function(X, ...) {
    normFactors = colSums(X)
    normFactors = normFactors / geom.mean(normFactors)
    return(sweep(X, 2, normFactors, "/"))
}

affy.loess <- function(X, refs.idx, verbose = FALSE, ...) {
    affy::normalize.loess(log2(X+1), subset = refs.idx, log.it = FALSE, verbose = verbose, ...)
}
