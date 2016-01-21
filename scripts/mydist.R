###############################################################################
### code chunk number : sampleClust--my sample distances function
###############################################################################
mydist = function (x, method="euclidean") {
        if (is.null(method) || missing(method)) {method <- "euclidean"}
        METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
                     "binary", "minkowski", "pearson", "spearman")
        method <- pmatch(method, METHODS)
        
        if (is.na(method)) stop("invalid distance method")
        if (any(method == seq(1,6))) {dist(t(x), method = METHODS[method])}
        else {as.dist(1-cor(x, method = METHODS[method]))}           
}