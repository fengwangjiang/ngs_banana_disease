###############################################################################
### code chunk number : MDS plots
###############################################################################
my.mds.plot <- function(rld, method="euclidean", dim=2, main=NULL, filename=NULL){
        library(MASS)
        if(dim !=2 && dim !=3) {
                stop(paste(dim, "can only be 2 or 3"))
        }
        meta <- as.data.frame(colData(rld))
        color <- c("red", "blue")
        legends <- levels(meta$condition)
        if (length(color) != length(legends))
                stop("length of color NOT equal to leves")
        jcolor <- data.frame(conditon=legends, color=I(color))
        
        label=colnames(rld)
        day <- paste(as.character(unique(meta$day)), collapse = "")
        if(is.null(main) || missing(main)){
                main <- paste("MDS ", dim, "D\n", method, " ", day, sep = "")
        }
        if(is.null(filename) || missing(filename)){
                filename <- paste0("mds_dim", dim, "_", method, day, ".pdf")
                filename <- file.path(DIR.FIGURE, filename)
        }
        
        # pdf(file = filename)
        
        # distance and isoMDS
        assay <- assay(rld)
        d <- mydist(assay, method)
        mds <- isoMDS(d, y=cmdscale(d, k=dim), k=dim)
        
        arglist_legend <- list(
                x = "topright",
                legend = as.character(jcolor$conditon),
                col = jcolor$color,
                pch = rep(16, length(legend)),
                bty = 'n',
                pt.cex = rep(2, length(legend))
        )
        if (dim == 2){
                my_df <- as.data.frame(mds$points)
                my_df <- dplyr::mutate(
                        my_df, name=rownames(my_df), x=V1, y=V2,
                        group=factor(ifelse(grepl("Cav", name), "Cav", "Bcl")),
                        cond=factor(ifelse(grepl("Ct", name), "Ct", "In")))
                        # group=as.factor(paste0(cell, cond)))
                # calculating the ellipses by df$group
                library(ellipse)
                # create an empty dataframe
                df_ell <- data.frame()
                for(g in levels(my_df$group)){
                        # create 100 points per variable around the mean of each group
                        df_ell <- rbind(df_ell, cbind(as.data.frame(with(my_df[my_df$group==g,], ellipse(cor(x, y), scale=c(sd(x),sd(y)), centre=c(mean(x),mean(y))))), group=g))
                }
                library(ggplot2)
                # create the ggplot with points colored by grouo
                g <- ggplot(data=my_df, aes(x=x, y=y,colour=group)) + 
                        # draw points
                        #geom_point(size=1.5, alpha=.6) +
                        geom_point(size=3, aes(shape=cond)) +
                        # draw ellipse lines
                        # geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=0.5, linetype=1) +
                        geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=0.5, linetype=0) +
                        # style as black and white theme
                        theme_bw()
                library(ggrepel)
                g_text <- g + geom_text_repel(data=my_df, aes(label=name))
                ggsave(filename = filename)
                
                # x = mds$points[, 1]
                # y = mds$points[, 2]
                # plot(x, y, col.axis="blue", type='p', cex=2, pch=16,
                #      xlim = 1.3*c(-max(abs(x)), max(abs(x))),
                #      ylim = 1.3*c(-max(abs(y)), max(abs(y))),
                #      main=main,
                #      col=jcolor$color[match(meta$condition, jcolor$conditon)],
                #      xlab="",ylab="",
                #      sub = sprintf("Stress = %.2f %%", mds$stress))
                # grid(col="lightblue")
                # wordcloud::textplot(x, y, words=label, cex=0.8, new = FALSE)
                # do.call(legend, arglist_legend)
        }
        if (dim == 3) {
                x = mds$points[,1]
                y = mds$points[,2]
                z = mds$points[,3]
                library(scatterplot3d)
                s3d<-scatterplot3d(x, y, z, 
                                   col.axis="blue", 
                                   col.grid="lightblue",
                                   color=jcolor$color[match(meta$condition, jcolor$conditon)],
                                   main=main,
                                   angle=60, type='h', cex.symbol=3, pch=16,
                                   highlight.3d=FALSE,
                                   xlab="", ylab="", zlab="",
                                   sub=sprintf("Stress = %.2f %%", mds$stress))
                coord.2d<-s3d$xyz.convert(x, y, z)
                wordcloud::textplot(coord.2d$x, coord.2d$y,
                                    words=label, cex=0.8, new = FALSE)
                do.call(legend, arglist_legend)
                tmp = dev.off()
        }
}