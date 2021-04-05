#This function takes in a vector of values and adjusts
#them to fit a particular range.

scale.between.vals <- function(V, target.min = 0, target.max = 1){

    r.min <- min(V, na.rm = TRUE)
    r.max <- max(V, na.rm = TRUE)

    new.vals <- ((target.max - target.min) * (V - r.min)/(r.max - r.min)) + target.min
    #plot(sort(new.vals))
    #plot(V, new.vals)
    return(new.vals)
}