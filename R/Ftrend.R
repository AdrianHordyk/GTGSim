#' Generate a trend in fishing mortality
#'
#' @param yr.st Numeric. First year - unfished
#' @param yr.end Numeric. Current year - F = `curF`
#' @param curF Numeric. Current apical fishing mortality
#' @param F.pat Character. F-pattern. Choices are: 'stable', 'inc' (increasing),
#' and 'dec' (decreasing)
#' @param Numeric. Fcv inter-annual CV of F
#' @param plot Logical. Should the F-trend be plotted?
#'
#' @return A numeric vector of fishing mortality
#' @export
#'
#' @examples
#' Ftrend(1950, 1990, 0.1, 'stable', 0.3, TRUE)
Ftrend <- function(yr.st, yr.end, curF,
                   F.pat= c('stable', 'inc', 'dec'),
                   Fcv=0.3, plot=TRUE) {
  F.pat <- match.arg(F.pat)
  yrs <- yr.st:yr.end
  yr.mid <- ceiling(mean(yrs))
  yr.ind <- which(yrs == yr.mid)
  nyrs <- length(yrs)
  Ferr <- exp(rnorm(nyrs, 0, Fcv))
  if (F.pat == "inc") {
    Ftrend <- seq(from=0, to=curF, length.out=nyrs) * Ferr
    Ftrend <- Ftrend/(Ftrend[length(Ftrend)]/curF)

  } else if (F.pat == "dec") {
    Ftrend <- c(seq(from=0, to=2*curF, length.out=yr.ind),
                seq(from=2*curF, to=curF, length.out=nyrs-yr.ind)) * Ferr
    Ftrend <- Ftrend/(Ftrend[length(Ftrend)]/curF)
  } else if (F.pat == "stable") {
    Ftrend <- c(seq(from=0, to=curF, length.out=yr.ind),
                seq(from=curF, to=curF, length.out=nyrs-yr.ind)) * Ferr
  }

  if (plot) {
    plot(yrs, Ftrend, type="l", bty="l", las=1, xlab="Years", ylab='apical Fishing mortality')
  }

  Ftrend
}
