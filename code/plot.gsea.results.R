#plots results from fgsea.


plot.gsea.results <- function(gsea.results, pname = "padj", p.thresh = 0.05, 
path.names = NULL, plot.label = NULL, label.cex = 0.8){
  if(nrow(gsea.results) == 0){
    plot.text("No significant results")
    return()
  }
  pcol <- which(colnames(gsea.results) == pname)
  if(is.null(path.names)){
    namecol <- which(colnames(gsea.results) == "pathway")
    path.names <- gsea.results[[namecol]]
  }
  
  pvals <- gsea.results[[pcol]]
  to.plot <- which(pvals <= p.thresh)

  if(length(to.plot) == 0){
    plot.text(paste("No signifiant results for", plot.label))
    return(NULL)
  }

  par(mar = c(4, 25,4, 4))
  result <- -log10(gsea.results[[pcol]][to.plot])
  result.names <- path.names[to.plot]
  result.order <- order(result)
  barplot(result[result.order], names = result.names[result.order], las = 2,
  xlab = paste("Negative log", pname), horiz = TRUE, cex.names = label.cex,
  main = plot.label)
}