
#' Make a Stock List with default values
#'
#' @return A list of Stock Values
#' @export
#'
#' @examples
#' MakeStock()
MakeStock <- function() {
  Stock <- list()
  Stock$M <- 0.2 # Annual natural mortality rate
  Stock$Linf <- 100
  Stock$LinfCV <- 0.15 # variability in Linf across GTGs
  Stock$K <- 0.1
  Stock$t0 <- -0.1 # must be < 0
  Stock$alpha <- 1E-5
  Stock$beta <- 3
  Stock$L50 <- 55
  Stock$L95 <- 65
  Stock$steepness <- 0.7
  Stock$sigmaR <- 0.6
  Stock$Linfsd <- 0.05  # inter-annual variability (log-normal CV)
  Stock$Ksd <- 0.05 # inter-annual variability (log-normal CV)
  Stock$t0sd <- 0.01 # inter-annual variability (log-normal CV)
  Stock$R0 <- 1000
  Stock
}

#' Make a Fleet List with default values
#'
#' @return A list of Fleet Values
#' @export
#'
#' @examples
#' MakeFleet()
MakeFleet <- function() {
  Fleet <- list()

  Fleet$yr.st <- 1950
  Fleet$yr.end <- 2017
  Fleet$F.pat <- "stable"

  Fleet$curF <- 0.1 # current apical F
  # selectivity-at-length - fishery
  Fleet$L5 <- 35
  Fleet$LFS <- 45
  Fleet$Vmaxlen <- 0.8 # selectivity at Linf
  Fleet

}

#' Make a Survey List with default values
#'
#' @return A list of Survey Values
#' @export
#'
#' @examples
#' MakeSurvey()
MakeSurvey <- function() {
  Survey <- list()
  # selectivity-at-length - survey
  Survey$sL5 <- 15
  Survey$sLFS <- 25
  Survey$sVmaxlen <- 1 # selectivity at Linf
  Survey$samp_ess <- 100 # annual effective sample size
  Survey$samp_n <- 300 # sample size
  Survey$ageSD <- 0.1
  Survey$lenSD <- 0.025
  Survey$wghtSD <- 0.05
  Survey
}


#' Simulate Data
#'
#' @param Stock A Stock List
#' @param Fleet A Fleet List
#' @param Survey A Survey List
#' @param Control A named list of control options
#' @param plot Logical. Show the plots?
#'
#' @return A data.frame of sampled age, length, weight, and maturity data
#' @export
#'
#' @examples
#' \dontrun{
#' MyData <- SimData()
#' }
SimData <- function(Stock=MakeStock(),
                    Fleet=MakeFleet(),
                    Survey=MakeSurvey(), Control=list(ngtg=101), plot=FALSE) {

  # unpack
  for (X in 1:length(Stock)) assign(names(Stock)[X], Stock[[X]])
  for (X in 1:length(Fleet)) assign(names(Fleet)[X], Fleet[[X]])
  for (X in 1:length(Survey)) assign(names(Survey)[X], Survey[[X]])
  for (X in 1:length(Control)) assign(names(Control)[X], Control[[X]])

  # Life-history
  maxage <- ceiling(-log(0.01)/M) # maxium age
  ages <- 0:maxage # in years

  # F pattern
  yrs <- yr.st:yr.end
  nyrs <- length(yrs)
  F.trend <- Ftrend(yr.st, yr.end, curF, F.pat, plot=plot)


  # selectivity-at-length - fishery
  sl <- (LFS - L5) /((-log(0.05,2))^0.5)
  sr <- (Linf - LFS) / ((-log(Vmaxlen,2))^0.5) # selectivity parameters are constant for all years

  # selectivity-at-length - survey
  sls <- (sLFS - sL5) /((-log(0.05,2))^0.5)
  srs <- (Linf - sLFS) / ((-log(sVmaxlen,2))^0.5) # selectivity parameters are constant for all years

  if (length(samp_ess) < nyrs) samp_ess <- rep(samp_ess, nyrs)[1:nyrs]
  if (length(samp_n) < nyrs) samp_n <- rep(samp_n, nyrs)[1:nyrs]


  # Mean length-at-age by year
  Linfs <- rlnorm(nyrs, log(Linf), Linfsd)
  Ks <- rlnorm(nyrs, log(K), Ksd)
  t0s <- -rlnorm(nyrs, log(-t0), t0sd)


  # Growth-type-groups
  distGTG <- seq(from=-3, to=3, length.out = ngtg)
  rdist <- dnorm(distGTG, 0, 1)/sum(dnorm(distGTG, 0, 1))
  ind <- as.matrix(expand.grid(1:nyrs, 1:(maxage+1),1:ngtg))
  Linfgtg <- Linfs[ind[,1]]+LinfCV*Linfs[ind[,1]]*distGTG[ind[,3]]

  L50gtg <- L95gtg <- array(NA, dim=c(nyrs, maxage+1, ngtg))
  L50gtg[ind] <- L50/Linf * Linfgtg  # assume constant L50/Linf but allow L50 to vary by year
  L95gtg[ind] <- L50gtg[ind] + (L95-L50)


  LAA <- array(NA, dim=c(nyrs, maxage+1, ngtg)) # Length-at-age by GTG and year
  LAA[ind] <- (Linfgtg * (1-exp(-Ks[ind[,1]] *(ages[ind[,2]] - t0s[ind[,1]]))))
  WAA <- alpha * LAA^beta  # Weight-at-age by GTG and year
  SAA <- array(dnormal(LAA, LFS, sl, sr), dim=c(nyrs, maxage+1, ngtg)) # selectivty-at-age by GTG and year
  sSAA <- array(dnormal(LAA, sLFS, sls, srs), dim=c(nyrs, maxage+1, ngtg)) # survey selectivty-at-age by GTG and year
  MAA <- 1/(1 + exp(-log(19) * ((LAA - L50gtg)/(L95gtg- L50gtg)))) # Maturity-at-age by GTG and year

  M_array <- array(M, dim=c(nyrs, maxage+1, ngtg))
  FAA <- array(NA, dim=c(nyrs, maxage+1, ngtg))
  FAA[ind] <- SAA * F.trend[ind[,1]] # fishing mortality at age by GTG and year
  ZAA <- FAA + M_array # Z-at-age by GTG and year

  # Unfished Year One
  Nunfished <- array(NA, dim=c(nyrs, maxage+1, ngtg))
  SB <- array(NA, dim=c(nyrs, maxage+1, ngtg))
  Nunfished[1,1,] <- rdist * R0 # distribute virgin recruitment
  Nunfished[1,2:(maxage+1),] <- matrix(Nunfished[1,1,], nrow=maxage, ncol=ngtg, byrow=TRUE) * exp(-apply(M_array[1,1:maxage,], 2, cumsum))

  SB[1,,] <- Nunfished[1,,] * WAA[1,,] * MAA[1,,]
  SB0 <- sum(SB[1,,])
  SBcurr <- Rec <- rep(NA, nyrs)
  SBcurr[1] <- SB0
  Rec[1] <- R0

  Nfished <- Nunfished
  recmu <- -0.5 * (sigmaR)^2
  recdevs <- exp(rnorm(nyrs, recmu, sigmaR))

  sample <- sample_det <- rep(list(), nyrs)
  message("Year 1 of ", nyrs)
  for (yr in 2:nyrs) {
    message("Year ", yr, " of ", nyrs)
    Rec[yr] <- BHSRR(SBcurr[yr-1], SB0, R0, steepness) # recruitment
    Nfished[yr,1,] <- Rec[yr] * recdevs[yr] * rdist

    Nfished[yr,2:(maxage+1),] <- Nfished[yr-1,1:(maxage),] * exp(-ZAA[yr-1,1:(maxage),])

    SB[yr,,] <- Nfished[yr,,] * WAA[yr,,] * MAA[yr,,]
    SBcurr[yr] <- sum(SB[yr,,])

    # Sample
    survAge <- Nfished[yr,,] * sSAA[yr,,]
    probsurvAge <- survAge/sum(survAge)
    sampind <- sample(1:length(probsurvAge), samp_ess[yr], replace=TRUE, prob=as.vector(probsurvAge))
    agearray <- matrix(ages, nrow=maxage+1, ncol=ngtg)
    SampAges_ess <- agearray[sampind]

    SampAges <- rep(SampAges_ess, samp_n[yr])[1:samp_n[yr]]
    SampAges <- SampAges + runif(length(SampAges), -.5, .5) # add sub-year variability to ages

    # generate sampled data
    linfGTG <- Linfs[yr] +LinfCV*Linfs[yr]*distGTG
    linfArray <- matrix(linfGTG, nrow=maxage+1, ncol=ngtg, byrow=TRUE)
    linfs <- rep(linfArray[sampind], samp_n[yr])[1:samp_n[yr]]
    Lens_det1 <- LAA[yr,,]
    Lens_det1 <- Lens_det1[sampind]
    Lens_det <- rep(Lens_det1, samp_n[yr])[1:samp_n[yr]]

    Lens <- linfs * (1-exp(-Ks[yr] *(SampAges - t0s[yr])))

    Wghts <- alpha * Lens ^ beta * exp(rnorm(samp_n[yr], 0, 0.3))
    L50s_det1 <- L50gtg[yr,,]
    L50s_det1 <- L50s_det1[sampind]
    L50s_det <- rep(L50s_det1, samp_n[yr])[1:samp_n[yr]]
    L95s_det1 <- L95gtg[yr,,]
    L95s_det1 <- L95s_det1[sampind]
    L95s_det <- rep(L95s_det1, samp_n[yr])[1:samp_n[yr]]
    L50s <- rep(L50s_det, samp_n[yr])[1:samp_n[yr]]
    L95s <- rep(L95s_det, samp_n[yr])[1:samp_n[yr]]
    Pmat <- 1/(1 + exp(-log(19) * ((Lens_det -L50s)/(L95s- L50s))))
    Maturity <- rbinom(samp_n[yr], 1, Pmat)

    # apply observation error
    SampAges_ob <- SampAges * exp(rnorm(samp_n[yr], -0.5 * (ageSD)^2, ageSD))
    Lens_ob <- Lens * exp(rnorm(samp_n[yr], -0.5 * (lenSD)^2, lenSD))
    Wght_ob <- Wghts * exp(rnorm(samp_n[yr], -0.5 * (wghtSD)^2, wghtSD))
    Mat_ob <- Maturity # assume no obs error

    if (samp_n[yr]>0) {
      sample[[yr]] <- data.frame(Year=yrs[yr], Age=SampAges_ob, Length=Lens_ob,
                                 Weight=Wght_ob, Maturity=Mat_ob)
      sample_det[[yr]] <- data.frame(Year=yrs[yr], Age=SampAges, Length=Lens,
                                 Weight=Wghts, Maturity=Maturity)
    }
  }
  SampDat <- do.call("rbind", sample)
  SampDat_det <- do.call("rbind", sample_det)
  out <- list()
  out$SampDat <- round(SampDat,2)
  out$SampDat_det <- round(SampDat_det,2)
  out
}

