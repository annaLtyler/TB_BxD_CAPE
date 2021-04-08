#' Plot interaction plot for traits and genetic markers
#' 
#' This internal function is called by plot.effects
#' to generate an interaction plot for two markers
#' relating to a trait. If marker2_vals is NULL,
#' the function instead shows the main effect for
#' marker1.
#' 
#' @param phenoV A vector of trait values 
#' @param marker1_vals A vector of genotype values 
#' for marker1
#' @param marker2_vals A vector of genotype values
#' for marker2.
#' @param pheno_name A string indicating the name of
#' the trait being plotted.
#' @param marker1_label A string indicating the name
#' of marker1
#' @param marker2_label A string indicating the name
#' of marker2
#' @param ymin A numeric value indicating the minimum 
#' y value for the plot. If NULL, it will be calculated
#' based on phenoV and error bars.
#' @param ymax A numeric value indicating the maximum
#' y value for the plot. If NULL, it will be calculated
#' based on phenoV and error bars.
#' @param error_bars A string indicating the type of error
#' bars to draw. Can be "sd" for standard deviation, "se"
#' for standard error, or "none".
#'
#' @return None
#' @importFrom stats interaction.plot
#' @keywords internal

plot_continuous_lines <- function(phenoV, marker1_vals, marker2_vals = NULL, 
pheno_name, marker1_label, marker2_label, prob = 0.95){
	
	oldPar <- par(no.readonly = TRUE)
	on.exit(oldPar)


	#show main effects if marker2 is NULL
	if(is.null(marker2_vals)){
		plot.with.model(marker1_vals, phenoV, xlab = marker1_label, 
        main = pheno_name, ylab = pheno_name)
	}
		
	#if we have values for two markers
	#set up the plot both ways
	if(!is.null(marker2_vals)){
			not_na <- which(!is.na(phenoV))
            df <- data.frame(cbind(marker1_vals[not_na], marker2_vals[not_na], phenoV[not_na]))
			colnames(df) <- c("Source", "Target", "Phenotype")
			fml <- as.formula("Phenotype ~ Source * Target")
			model <- lm(fml, data = df)
			fig1 <- interact_plot(model, Source, Target, interval = TRUE, int.width = prob)
			fig2 <- interact_plot(model, Target, Source, interval = TRUE, int.width = prob)
			fig.list <- list(fig1, fig2)
			return(fig.list)
	}#end case for if there are two markers
		
}#end function
		