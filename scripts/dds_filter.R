dds_filter <- function(dds, low=10, high=1e4) {
        cnts <- counts(dds)
        rmin <- apply(cnts, 1, min)
        rmax <- apply(cnts, 1, max)
        if (is.null(low) || missing(low)) {
                low = max(10, quantile(cnts, probs=0.10))
        }
        if (is.null(high) || missing(high)) {
                high <- min(1e4, quantile(cnts, probs=0.99))
        }
        chosen <- rmin > 10 & rmax < high
        dds_chosen <- dds[chosen, ]
        return(dds_chosen)
}

# dds2 <- dds_filter(dds)
# length(dds2) / length(dds) # 0.5644346