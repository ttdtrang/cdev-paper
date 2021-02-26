#' Geometric mean
geom.mean <- function(x) {
    log(x) %>%
        mean(na.rm = TRUE) %>%
        exp() %>%
        return()
}

#' Relative log expression
relative_log_expression <- function(x, log.base = 2, pseudocount = 1) {
    logx = log(x+pseudocount, base = log.base)
    med = apply(logx, MARGIN = 1, median)
    sweep(logx, 1, STATS = med)
}

cv <- function(x) {
    return(sd(x)/mean(x))
}

norm.minmax <- function(x, na.rm = TRUE) {
    return((x-min(x, na.rm = TRUE))/(max(x, na.rm = na.rm) - min(x, na.rm = na.rm)))
}
