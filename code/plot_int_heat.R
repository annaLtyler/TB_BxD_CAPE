#' Plot phenotypic effects for two markers as a heat map
#' 
#' This internal function is called by 
#' \code{\link{plot_effects}} to generate a 
#' heat map showing the effects of genotype on
#' phenotype. This function fits linear models
#' to the markers and traits. It then uses
#' these models to predict trait values at 
#' different genotype combinations in a 2D 
#' grid. It plots these predicted values in
#' a heat map.
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
#' @param bins1 The number of bins for marker1 over 
#' which to predict values of the trait. This can also 
#' be a vector specifying those bins.
#' @param bins2 The number of bins for marker2 over 
#' which to predict values of the trait. This can also 
#' be a vector specifying those bins.
#'
#' @return None
#' 
#' @importFrom stats as.formula lm predict
#' @keywords internal

plot_int_heat <- function(phenoV, marker1_vals, marker2_vals, pheno_name = NULL, 
    marker1_label = NULL, marker2_label = NULL, bins1 = 50, bins2 = 50){

	has_dash1 <- grep("-", marker1_label)
	if(length(has_dash1) > 0){
		marker1_label <- gsub("-", "_", marker1_label)
	}
	has_dash2 <- grep("-", marker2_label)
	if(length(has_dash2) > 0){
		marker2_label <- gsub("-", "_", marker2_label)
	}
	

    if(length(bins1) == 1){
	    min_marker1 <- min(signif(marker1_vals, 2), na.rm = TRUE)
		max_marker1 <- max(signif(marker1_vals, 2), na.rm = TRUE)
        marker_grid1 <- segment_region(min_marker1, max_marker1, bins1, alignment = "ends")
    }else{
    	marker_grid1 <- bins1
    }
    
    if(length(bins2) == 1){
		min_marker2 <- min(signif(marker2_vals, 2), na.rm = TRUE)
		max_marker2 <- max(signif(marker2_vals, 2), na.rm = TRUE)
        marker_grid2 <- segment_region(min_marker2, max_marker2, bins2, alignment = "ends")
    }else{
       	marker_grid2 <- bins2
   	}
        
    marker1_bins <- bin_vector(signif(marker1_vals, 2), bins = marker_grid1)
    marker2_bins <- bin_vector(signif(marker2_vals, 2), bins = marker_grid2)
    #plot(marker1_bins, marker2_bins)

	test_df <- data.frame(phenoV, marker1_bins, marker2_bins)
	colnames(test_df) <- c(pheno_name, marker1_label, marker2_label)
	int.fmla <- paste(pheno_name, "~", marker1_label, "*", marker2_label)
    int.model <- lm(as.formula(int.fmla), data = test_df)
	add.fmla <- paste(pheno_name, "~", marker1_label, "+", marker2_label)
	add.model <- lm(as.formula(add.fmla), data = test_df)
	#summary(add.model)
	
	predict_grid <- cbind(rep(marker_grid1, length(marker_grid2)), rep(marker_grid2, 
	each = length(marker_grid1)))
	colnames(predict_grid) <- c(marker1_label, marker2_label)
	predict_df <- data.frame(predict_grid)
	
	pred_data_add <- predict(add.model, newdata = predict_df, se.fit = TRUE)
	pred_mat_add <- matrix(pred_data_add$fit, nrow = length(marker_grid2), 
	ncol = length(marker_grid1), byrow = TRUE)
	rownames(pred_mat_add) <- signif(marker_grid2, 2)
	colnames(pred_mat_add) <- signif(marker_grid1, 2)
	#pred_cent_add <- pred_mat_add - mean(pred_mat_add)
	pred_add_se <- matrix(pred_data_add$se.fit, nrow = length(marker_grid2), 
	ncol = length(marker_grid1), byrow = TRUE)
	
	pred_data_int <- predict(int.model, newdata = predict_df, se.fit = TRUE)
	pred_mat_int <- matrix(pred_data_int$fit, nrow = length(marker_grid2), 
	ncol = length(marker_grid1), byrow = TRUE)
	rownames(pred_mat_int) <- signif(marker_grid2, 2)
	colnames(pred_mat_int) <- signif(marker_grid1, 2)
	#pred_cent_int <- pred_mat_int - mean(pred_mat_int)
	pred_int_se <- matrix(pred_data_int$se.fit, nrow = length(marker_grid2), 
	ncol = length(marker_grid1), byrow = TRUE)

	#pred_cent_diff <- pred_cent_int - pred_cent_add
	pred_diff <- pred_mat_int - pred_mat_add

	all.min <- min(sapply(list(pred_mat_add-pred_add_se, pred_mat_int-pred_int_se, 
	pred_diff), min))
	all.max <- max(sapply(list(pred_mat_add+pred_add_se, pred_mat_int+pred_int_se, 
	pred_diff), max))

	par(mfrow = c(1,3))

	image_with_text(pred_mat_add[nrow(pred_mat_add):1,], col_text_rotation = 0, 
	use_pheatmap_colors = TRUE, show_text = FALSE, ylab = marker2_label, 
	xlab = marker1_label, global_color_scale = TRUE, global_min = all.min, 
	global_max = all.max,main = "Additive")
	#boxplot(pred_cent_add, main = "Marginal Effects of Source")
	#boxplot(t(pred_cent_add), main = "Marginal Effects of Target")

	image_with_text(pred_mat_int[nrow(pred_mat_int):1,], col_text_rotation = 0, 
	use_pheatmap_colors = TRUE, 
	show_text = FALSE, ylab = marker2_label, xlab = marker1_label,
	global_color_scale = TRUE, global_min = all.min, global_max = all.max, 
	main = "Interaction")
	#boxplot(pred_cent_int, main = "Marginal Effects of Source")
	#boxplot(t(pred_cent_int), main = "Marginal Effects of Target")

	image_with_text(pred_diff[nrow(pred_diff):1,], col_text_rotation = 0, 
	use_pheatmap_colors = TRUE, show_text = FALSE, ylab = marker2_label, 
	xlab = marker1_label, global_color_scale = TRUE, global_min = all.min, 
	global_max = all.max, main = "Difference")

	#image_with_text(pred_int_se, col_text_rotation = 0, use_pheatmap_colors = TRUE, 
	#show_text = FALSE, ylab = marker2_label, xlab = marker1_label,
	#main = "SE")
	
	colorscale_add <- list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)"))
	colorscale_int <- list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))

	fig <- plot_ly(showscale = FALSE)
	fig <- add_surface(fig, z = ~pred_mat_add, opacity = 0.90, type = "surface", 
	colorscale = colorscale_add, showlegend = TRUE, name = "Additive")
	fig <- add_surface(fig, z = ~pred_mat_add+pred_add_se, opacity = 0.50, 
	type = "surface", colorscale = colorscale_add)
	fig <- add_surface(fig, z = ~pred_mat_add-pred_add_se, opacity = 0.50, 
	type = "surface", colorscale = colorscale_add)

	fig <- add_surface(fig, z = ~pred_mat_int, opacity = 0.90, type = "surface",
	colorscale = colorscale_int, showlegend = TRUE, name = "Interactive")
	fig <- add_surface(fig, z = ~pred_mat_int-pred_int_se, opacity = 0.50, 
	type = "surface", colorscale = colorscale_int)
	fig <- add_surface(fig, z = ~pred_mat_int+pred_int_se, opacity = 0.50, 
	type = "surface", colorscale = colorscale_int)
	#fig
	return(fig)
}