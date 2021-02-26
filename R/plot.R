#' plot a matrix into heatmap
#'
#' @import ggplot2
plot.heatmap <- function(X) {
    mm = reshape2::melt(X)
    p = ggplot(mm) +
        geom_tile(aes(x=Var1,y=Var2,fill=value)) +
        scale_fill_gradient2(low='red', mid='white', high='steelblue') +
        theme(axis.text.x = element_text(angle=30,hjust=1), axis.title.x = element_blank(), axis.title.y = element_blank()) +
        coord_fixed()
    return(p)
}

#' plot an expression matrix as parallel coordinate plot
#'
#' @param X a matrix with variables in columns and observations in rows
#' @param alpha a transparency value for lines, ranging from 0 to 1
#' @param ylabel Label of the y axis
#' @param xlabel Label of the x axis
#' @param line.color A scalar or vector specifying how each line should be colored
#' @import ggplot2, reshape2
plot.pcp <- function(X, alpha = 0.3, ylabel='value', xlabel='Var2', line.color = 'darkgreen') {
    meltedExpr = reshape2::melt(X)

    if (length(line.color) > 1) {
        df.color = data.frame('Var1' = rownames(X),
                              'line.color' = line.color)
        meltedExpr = merge(x = meltedExpr, y = df.color, by = 'Var1', all.x = TRUE, all.y=FALSE)
        the_line_geom = geom_line(aes(color=line.color),alpha= alpha)
    } else {
        the_line_geom = geom_line(alpha= alpha,color = line.color)
    }


    p = ggplot(meltedExpr, aes(x=Var2,y=value,group=Var1)) +
        the_line_geom +
        scale_y_continuous(trans = 'log10') +
        theme(axis.text.x = element_text(angle=60,hjust=1),legend.position='bottom') +
        ylab(ylabel) +
        xlab(xlabel)
    return(p)
}

#' plot a raster image of a matrix
#'
plot.matrix = function(m, col = colorRampPalette(colors = c("blue", "white", "red"))(n = 999), asp=1, ...) {
    # burd = colorRampPalette(colors = c("blue", "white", "red"))(n = 999)
    # blues = colorRampPalette(colors = c('white', 'blue'))(n = 255)
    m %>%
        apply(MARGIN = 2, rev) %>%
        t() %>%
        image(useRaster = TRUE, axes = FALSE, col = col, asp = asp, ...)
}
