#This function takes an x and y and plots them with a fitted linear model


plot.with.model <- function(x, y, xlim = NULL, ylim = NULL, col = "black", 
pch = 16, main, xlab, ylab, report = c("lm", "cor.test"), cex = 1,
plot.results = TRUE, add = FALSE, plot.type = "p"){
	
	if(missing(main)){main = ""}
	if(missing(xlab)){xlab = deparse(substitute(x))}
	if(missing(ylab)){ylab = deparse(substitute(y))}

	report <- report[1]
	
	model <- lm(y~x)
	
	if(report == "lm"){
		f <- summary(model)$fstatistic
	    p <- pf(f[1],f[2],f[3],lower.tail=F)
		r2 <- signif(summary(model)$r.squared, 2)
		rlab <- "R2"
		}else{
		test <- cor.test(x,y)	
		p <- signif(test$p.value, 2)
		r2 <- signif(test$estimate, 2)
		rlab <- "r"
		}

	if(plot.results){
	new.title <- paste0(main, "\n", rlab, " = ", r2, ", p = ", signif(p, 2))
	if(!add){
		plot(x,y, xlim = xlim, ylim = ylim, col = col, pch = pch, main = new.title, 
		xlab = xlab, ylab = ylab, type = plot.type)
	}else{
		points(x,y, col = col, pch = pch)
	}
	abline(model, col = "#a6bddb", lwd = 3)
	}
	
	invisible(c("r" = r2, "p" = p))
	
}