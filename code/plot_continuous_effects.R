plot_continuous_effects  <- function(data_obj, geno_obj, marker1, marker2 = NULL, 
	pheno_type = "normalized", plot_type = c("l", "h"), sig.dig = 3,
	prob = 0.95, covar = NULL, marker1_label = NULL, marker2_label = NULL, 
    bins_marker1 = 50, bins_marker2 = 50, separate_windows = FALSE){


    plot_type = plot_type[1]
		
	#==========================================
	# get traits for plotting
	#==========================================
	covar_info <- get_covar(data_obj)		
	
	if(is.null(covar)){
		covar_names <- covar_info$covar_names
	}else{
		if(covar == "none"){
			covar_names = NULL
		}else{
			covar_names = covar
		}
	}
	
    pheno <- get_pheno(data_obj, pheno_type, covar_names)
	#==========================================


	#==========================================
	# get marker genotypes or covariate values 
	# for plotting
	#==========================================
	marker_vals <- get_marker_covar(data_obj, geno_obj, c(marker1, marker2))
	marker_vals <- apply(marker_vals, 2, function(x) signif(x, sig.dig))


   	#==========================================
	# line up genotypes and phenotypes
	#==========================================
    common_ind <- intersect(rownames(marker_vals), rownames(pheno))
    ind_pheno_locale <- match(common_ind, rownames(pheno))
    ind_geno_locale <- match(common_ind, rownames(marker_vals))
    
    pheno_to_plot <- pheno[ind_pheno_locale,,drop=FALSE]
    geno_to_plot <- marker_vals[ind_geno_locale,,drop=FALSE]
 	#==========================================


	#============================================================
	# assign the marker names
	#============================================================
	if(is.null(marker1_label)){marker1_label = marker1}
	if(is.null(marker2_label)){marker2_label = marker2}
	marker_names <- c(marker1_label, marker2_label)
	#============================================================

	#============================================================
	# get layout matrix for heat maps
	# the interaction plots use ggplot and cannot be put on 
	# the same device
	#============================================================
	layout.mat <- get.layout.mat(ncol(pheno))
	if(plot_type == "h"){
		layout(layout.mat)
	}
	#============================================================	
	fig.list <- vector(mode = "list", length = ncol(pheno))
	for(ph in 1:ncol(pheno_to_plot)){
		if(separate_windows){
			dev.new()
		}
		phenoV = pheno_to_plot[,ph]
		pheno_name = colnames(pheno)[ph]
		marker1_vals <- geno_to_plot[,1]
		if(!is.null(marker2)){
			marker2_vals <- geno_to_plot[,2]
		}else{
			marker2_vals <- NULL
		}
		
		if(plot_type == "h"){
			if(is.null(marker2_vals)){stop("Two markers are required for the heat map.")}
		  	fig.list[[ph]] <- plot_int_heat(phenoV, marker1_vals, marker2_vals, pheno_name,
			marker1_label, marker2_label, bins1 = bins_marker1, bins2 = bins_marker2)
		}
		if(plot_type == "l"){
		  fig.list[[ph]] <- plot_continuous_lines(phenoV, marker1_vals, marker2_vals, pheno_name, 
			marker1_label, marker2_label, prob)
		}	
		
	} #end looping through phenotypes


	return(fig.list)
	

    }