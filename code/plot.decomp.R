#This function mean centers a matrix,
#takes the svd and then plots the specified
#principal components against each other.
#cols is a vector of color values for the points
#color.by is a numeric matrix of factors whose
#relationship to the principal components can 
#be examined. 
#color.by supercedes cols if specified.

plot.decomp <- function(mat, pc = 2, mean.center = TRUE, 
scale.col = FALSE, cols = rep("black", nrow(mat)), color.by = NULL, 
pch = NULL, cex = 1, main = NULL, plot.results = TRUE){

	na.rows <- which(apply(mat, 1, function(x) all(!is.na(x))))
	if(length(na.rows) < nrow(mat)){
		warning("Removing rows with missing values")
	}

	if(mean.center){
		cent <- apply(mat[na.rows,], 2, function(x) mean(x, na.rm = TRUE) - x)
	}else{
		cent <- mat
	}

	if(scale.col){
		scaled.mat <- scale(cent)
	}else{
		scaled.mat <- cent
	}
	
	decomp <- svd(scaled.mat, nu = pc, nv = pc)
	pc.mat <- decomp$u
	var.exp <- decomp$d^2/sum(decomp$d^2)
	var.text <- signif(var.exp[1:ncol(pc.mat)]*100, 2)	
	colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat), " (", var.text, "%)")
	decomp$rows.used <- na.rows
	
	if(plot.results){
		if(is.null(pch)){pch = 16}

		if(is.null(color.by)){
			num.plots = 1
			plot.label <- main
		}else{
			num.plots <- ncol(color.by)
			if(is.null(colnames(color.by))){
				colnames(color.by) <- paste(main, "Factor", 1:ncol(color.by))
			}
			plot.label <- colnames(color.by)
		}
			
		for(i in 1:num.plots){
			if(!is.null(color.by)){
				if(is.numeric(color.by[na.rows,i])){
					cols <- colors.from.values(color.by[na.rows,i])
				}else{
					cols <- as.numeric(as.factor(color.by[na.rows,i]))
				}
			}
			if(pc == 2){
				plot(pc.mat[,1], pc.mat[,2], xlab = colnames(pc.mat)[1], 
				ylab = colnames(pc.mat)[2], 
				main = plot.label[i], col = cols, pch = pch, cex = cex)
			}else{
				pairs(pc.mat[,1:pc], col = cols, pch = pch, cex = cex,
				main = plot.label[i])
				}
		}
	}
	
	invisible(decomp)

}
