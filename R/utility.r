## facility functions ####
Gmat <- function(nvec){
  D <- length(nvec)
  arealab <- rep(1:D, times = nvec)
  unname(model.matrix(~as.factor(arealab)-1))
}
