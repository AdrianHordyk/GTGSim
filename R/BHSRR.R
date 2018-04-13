BHSRR <- function(SBcurr, SB0, R0, steepness) {
  (4 * R0 * steepness * SBcurr)/(SB0/R0 * R0 * (1-steepness) + (5*steepness-1)*SBcurr)
}
