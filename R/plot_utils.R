
#' Generate the defult colors
#'
#' @param  n integer, number of colors
#'
#' @return n colors
#' @export
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#' UMAP plot
#'
#' @param X n-by-p expression matrix
#' @param labels vector of sample labels
#' @param pca umap parameter, dimension of pca
#' @param n_components umap parameter, dimension of low dimensional embedding space (default 2)
#' @param n_neighbors umap parameter, number of neighbors for nearest neighbor graph
#' @param min_dist umap parameter, minimum distance of point to its nearest neighbor in the
#' embedding space.
#' @param point.size numeric scalar, point size in the plot
#' @param alpha numeric, transparency of the points in the plot
#' @param title string, title of the plot
#' @param legend.name string, legend name
#' @param cols vector of colors, length should the same as cardinality of labels
#' @param emb embedding of the UMAP if provided
#' @param seed random seed
#'
#' @return A list of \describe{
#'  \item{p}{umap plot}
#'  \item{emb}{umap embedding matrix}
#' }
#'
#' @import ggplot2
#' @importFrom uwot umap
#' @export
plot_umap <- function(X = NULL, labels = NULL, pca = 50, n_components = 2, n_neighbors = 30, 
    min_dist = 0.1, point.size = 0.3, alpha = 1, title = NULL, legend.name = "labels", 
    cols = NULL, emb = NULL, seed = 0) {
    library(ggplot2)
    
    if (is.null(X) & is.null(emb)) {
        stop("data not provided!")
    }
    
    set.seed(seed)
    
    if (is.null(emb)) {
        if (!is.null(pca)) {
            if (pca > ncol(X)/2) {
                pca = NULL
            }
        }
        emb = uwot::umap(X, n_neighbors = n_neighbors, n_components = n_components, 
            min_dist = min_dist, pca = pca)
    }
    
    df = data.frame(umap1 = emb[, 1], umap2 = emb[, 2], labels = if (!is.null(labels)) 
        labels else rep(0, nrow(X)))
    p = ggplot2::ggplot(df, aes(x = umap1, y = umap2)) + geom_point(col = "black", 
        size = point.size, stroke = 0, shape = 16, alpha = alpha) + labs(x = "UMAP_1", 
        y = "UMAP_2", title = title) + theme_light() + theme(plot.title = element_text(hjust = 0.5), 
        axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
    
    if (!is.null(labels)) {
        if (is.null(legend.name)) {
            legend.name = "labels"
        }
        
        if (is.null(cols)) {
            cols = gg_color_hue(length(unique(labels)))
        }
        
        p = p + geom_point(aes(colour = labels), size = point.size, stroke = 0, shape = 16, 
            alpha = alpha) + scale_color_manual(values = cols) + guides(col = guide_legend(ncol = 1, 
            title = legend.name, override.aes = list(size = 5)))
    }
    list(p = p, emb = emb)
}
