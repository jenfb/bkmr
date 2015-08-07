##http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

reorder_cormat <- function(cormat) {
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <- cormat[hc$order, hc$order]
    return(cormat)
}

# Get lower triangle of the correlation matrix
get_lower_tri <-function(cormat) {
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
}

#' Plot of correlation matrix
#'
#' @export
#' @import ggplot2
#' @param vals
#' @param digits If the correlation values are to be printed in the plot (\code{vals == TRUE}), how many digits should be displayed
PlotCorr <- function(Z, reorder = TRUE, print.vals = FALSE, plot = TRUE, digits = 2) {
    cormat <- cor(Z)
    if(print.vals) cormat <- round(cormat, digits)

    # Reorder the correlation matrix
    if(reorder) cormat <- reorder_cormat(cormat)
    upper_tri <- get_upper_tri(cormat)

    # Melt the correlation matrix
    melted_cormat <- reshape2::melt(upper_tri)
    melted_cormat <- na.omit(melted_cormat)

    # Create a ggheatmap
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                             midpoint = 0, limit = c(-1,1), name="Pearson\nCorrelation") +
        theme_minimal() + # minimal theme
        theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                         size = 12, hjust = 1),
              axis.text.y = element_text(size = 12)) +
        coord_fixed() +
        labs(x = "", y = "")
    if(print.vals) {
        ggheatmap <- ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 3)
    }

    # Print the heatmap
    if(plot) print(ggheatmap)
    invisible(ggheatmap)
}



