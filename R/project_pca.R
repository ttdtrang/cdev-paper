#' project.pca
#' 
#' Project a given set of data points on 
#' selected principal components
#' @param X data matrix, points in rows, dimensions in columns. In case of gene expression matrix, it would likely mean samples x genes.
#' @param selected.pc all PCs are selected by default
#' @param center whether or not to center the data to zero mean before SVD
#' @param scale whether or not to scale the data to unit variance before SVD
#' @param sv_scale whether or not to scale the projected data points by their corresponding singular values
#' @export
project.pca <- function(X, selected.pc = c(), center = TRUE, scale = FALSE, sv_scale = FALSE) {
    s = X %>%
        scale(center = center, scale = scale) %>%
        svd()
    if (length(selected.pc) == 0) selected.pc = 1:length(s$d)
    projected = s$u[,selected.pc]
    if (sv_scale) {
        projected = projected %*% diag(s$d[selected.pc])
    }
    data.frame(projected, row.names = rownames(X), check.names = FALSE) %>%
        set_names(paste0('PC', selected.pc)) %>%
        return()
}
