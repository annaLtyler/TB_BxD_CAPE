#' Exponential color function
#' 
#' This function maps colors onto an exponential
#' function for use in colors_from_values.
#' 
#' @param x.min The minimum value assigned to a color
#' @param x.max The maximum value assigned to a color
#' @param steepness The steepness of the exponential function
#' @param num.cols The number of colors to generate in the ramp
#' 
#' @return Numeric vector mapping x.min:x.max onto an
#' exponential curve 

exp.color.fun <-
function(x.min, x.max, steepness = 1, num.cols = 256){
	
	x <- seq(0,1,length.out = num.cols)
	y <- exp(steepness*x)-1
	
	#map the 0-1 interval onto the interval from x.min to x.max
	x.range <- x.max - x.min
	scaled.y <- (y*x.range)/max(y)
	

	shifted.y <- scaled.y + x.min
	
	return(shifted.y)
	
	
}
