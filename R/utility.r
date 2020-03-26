## facility functions ####
expo <- function(x) ifelse(x!=0, exp(x), 0)
logo <- function(x) ifelse(x>0, log(x), 0)
Gmat <- function(nvec){
  D <- length(nvec)
  arealab <- rep(1:D, times = nvec)
  unname(model.matrix(~as.factor(arealab)-1))
}
