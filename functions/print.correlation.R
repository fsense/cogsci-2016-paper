#
# Florian Sense [f.sense@rug.nl]
#
# Function to print details about a correlation to the console
#
print.correlation <- function(x, y, x.label="variable A", y.label="variable B") {
  if(length(x) != length(y)) stop("x and y must be equally long!")
  tmp <- cbind(x, y)
  tmp <- na.omit(tmp) # make sure there are no NAs

  corr <- cor.test(tmp[, 1], tmp[, 2])

  # Function from Wetzels & Wagenmakers (2012):
  cor.bf=function(r,n){
    int=function(r,n,g){
      exp(
        ((n-2)/2)*log(1+g)+
          (-(n-1)/2)*log(1+(1-r^2)*g)+
          (-3/2)*log(g)+
          -n/(2*g))
    }
    bf10=sqrt((n/2))/gamma(1/2)*
      integrate(int,lower=0,upper=Inf,r=r,n=n)$value
    return(bf10)
  }

  bf <- cor.bf(corr$estimate, nrow(tmp))

  cat(paste("cor(", x.label, ", ", y.label, ") = ", round(corr$estimate, 3), " with t(", corr$parameter, ") = ", round(corr$statistic, 2), " and p >= ", round(corr$p.value, 4), ".\n    The BF is ", ifelse(bf < 1, paste(round(1/bf, 2), "(in favor of null)."), paste(round(bf, 2), "(in favor of alt.).")) , sep=""))
}