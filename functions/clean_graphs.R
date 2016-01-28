###
### Florian Sense [f.sense@rug.nl]
###
### A collection of simple, clean graphs.
###
### Most of these are based on http://shinyapps.org/apps/RGraphCompendium/index.php#including-a-density-estimator
###

# colors for plotting
# Taken from: http://flatuicolors.com/#
carrot <- "#e67e22"
river <- "#3498db"
silver <- "#bdc3c7"
col.correct <- "#27ae60"
col.incorrect <- "#c0392b"
spancol <- c("#27ae60", "#2980b9", "#8e44ad") # for the individual complex span tasks

### Histogram:
clean.hist <- function(x, xlab="", ylab="", xlim=NULL, main=NULL, color=river, plot.density=TRUE, plot.rug=FALSE, rug.jitter=.6, shrink.mar=TRUE, breaks="Sturges", return.details = FALSE) {

  par(las=1)
  if(is.null(main) & shrink.mar) {
    # remove the top margin from the plotting area:
    op <- par()$mar
    op[3] <- .1 # make top margin small
    par(mar=op)
  }

  if(is.null(xlim)) {
    x <- x[!is.na(x)]
    if(max(x) < 10) {
      # probably not perfect but should work for most smaller numbers:
      xlim <- c(floor(min(x, na.rm = TRUE)), ceiling(max(x, na.rm = TRUE)))
    } else {
      xlim <- range(x) * c(.9, 1.1) # this kind of works most of the time
    }
  }

  # the plot:
  y <- hist(x, xlab=xlab, ylab=ylab, xlim=xlim, main=main, col=color, freq=ifelse(plot.density, FALSE, TRUE), breaks=breaks, axes=FALSE, border=NA)
  axis(1)

  if(plot.rug) {
    rug(jitter(x, rug.jitter))
  }

  if(plot.density) {
    lines(density(x[!is.na(x)]), lwd=2)
  }

  if(return.details) return(y)
}


### Scatter plot (with regression line)
clean.plot <- function(x, y, plot.regr = TRUE, main=NULL, xlab="", ylab="", ylim=NULL, xlim=NULL, color=river, asp=NULL, shrink.mar=TRUE) {

  # remove NAs if necessary:
  y <- y[!is.na(y)]
  x <- x[!is.na(x)]

  par(las=1, bty='n')
  if(is.null(main) & shrink.mar) {
    # remove the top margin from the plotting area:
    op <- par()$mar
    op[3] <- .2 # make top margin small
    par(mar=op)
  }

  plot(x, y, col=color, pch=20, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, main=main, asp=asp)

  if(plot.regr) {
    library(plotrix) # for the ablineclip() function
    regr <- lm(y ~ x)
    ablineclip(regr, lwd=1.5, x1=min(x), x2=max(x))
  }
}


### Overview plot (grid with scatter plots, densities, and correlations)
clean.overview <- function(dat, variables=NULL, colors=NULL, print.N=TRUE, print.bf=TRUE) {
  #
  # dat: a data.frame in which each column corresponds to one of the variables
  #
  if(print.bf) {
    cor.bf=function(r,n){
      # Function taken from Wetzels & Wagenmakers (2012);
      # code obtained from Ruud Wetzel's website
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
  }

  # Determine some defaults:
  if(is.null(variables)) {
    if(is.null(colnames(dat))) {
      variables <- LETTERS[1:ncol(dat)]
    } else {
      variables <- colnames(dat)
    }
  }

  if(is.null(colors)) {
    colors <- 1:ncol(dat) # this is ugly but will work.
  }

  if(ncol(dat) != length(colors) | ncol(dat) != length(variables)) {
    stop("Input lengths of either `colors` or `variables` doesn't match number of columns in `dat`!")
  }

  l <- length(variables)

  par(mfrow= c(l,l), cex.axis= 1.3, mar= c(3, 4, 2, 1.5) + 0.1, oma= c(0, 2.2, 2, 0))

  for(y in 1:length(variables)) {
    for(x in 1:length(variables)) {

      # Make sure all cases are complete:
      tmp <- as.data.frame(na.omit(cbind(dat[, x], dat[, y])))
      colnames(tmp) <- c("x", "y")

      # plot distr. + density on the diagonal
      if(y == x) {
        xlim <- range(tmp$y)
        out <- clean.hist(tmp$y, col=colors[y], xlim=xlim, shrink.mar=FALSE, plot.rug = TRUE, plot.density=FALSE, return.details=TRUE)

        if(print.N) {
          text(mean(range(tmp$y)), max(out$counts)/8, paste("N =", length(tmp$y)), cex=1.3)
        }
      } else if(x > y) {
        clean.plot(tmp$x, tmp$y, col="grey") #, xlim=c(0, max(dat[, x])))
      } else {
        plot(1, 1, type='n', axes=FALSE, xlab="", ylab="")
        corr <- cor(tmp$x, tmp$y)
        text(1, 1.1, paste("r =", round( corr, 2 )), cex=2.1)
        text(1, .9 , paste("p >=", round( cor.test(tmp$x, tmp$y)$p.value, 3 )), cex=1.5)
        if(print.bf) {
          bf <- cor.bf(corr, nrow(tmp))
          if(bf > 1000000 | bf < 1/1000000) {
            # Once it gets this large, the precise number doesn't really matter any more (in most cases):
            txt <- "BF > 1M"
          } else if(bf > 10000 | bf < (1/10000)) {
            # round down to nearest 1,000:
            txt <- paste("BF > ", floor(ifelse(bf > 1, bf, 1/bf)/1000), "k", sep="")
          } else {
            # This is the range in which we want to display the actual numbers:
            if(bf > 1) {
              txt <- bquote(BF['H1'] ~ "=" ~ .(round(bf, ifelse(bf > 100, 0, 1))))
            } else {
              txt <- bquote(BF['H0'] ~ "=" ~ .(round(1/bf, ifelse((1/bf) > 100, 0, 1))))
            }
          }

          text(1, .75 , txt, cex=1.5)
        }
      }
    }
  }

  textpos <- seq(1/(l*2), (l*2-1)/(l*2), 2/(l*2))
  for(t in seq_along(textpos)){
    mtext(text = variables[t], side = 3, outer = TRUE, at= textpos[t], cex=1.2, line= .3)
    mtext(text = variables[t], side = 2, outer = TRUE, at= rev(textpos)[t], cex=1.2, line= -0.4, las=0)
  }
}