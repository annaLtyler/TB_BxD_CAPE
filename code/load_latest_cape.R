#This function loads a local version of cape code.

load_latest_cape <- function(cape.dir){

    source(here("code", "load_libraries.R"))
    needed.libraries <- c("here", "qtl2", "abind", "Matrix", "MASS", "regress", "igraph",
        "RColorBrewer", "R6", "yaml", "tools", "caTools", "propagate", "igraph",
        "qtl", "regress", "evd")
    load_libraries(needed.libraries)

    cape.fun <- list.files(cape.dir, pattern = ".R", full.names = TRUE)
    for(i in 1:length(cape.fun)){source(cape.fun[i])}

}