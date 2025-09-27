library(ggplot2)
library(gridExtra)
library(latex2exp)
nVec = 16:95
numNVals = length(nVec)
stdDevVec = c(3,3,4.5,4.5)
muVec = c(70,52,70,52)
for (iPlot in 1:4) {
  mu = muVec[iPlot]
  s = stdDevVec[iPlot]
  nTermsVec = c(4,7,10)
  pVec = c(4,7,10)
  t = 2^3 + 3*2 + 1
  powerMat = matrix(0,numNVals,3)
  iCol = 0
  # for (p in pVec) {
  p = 10
  fVec = c(.025,.05,0.1)
  plotPointsDF = matrix(0,5,4)
  for (f in fVec) {
    iCol = iCol + 1
    numDF = t - p
    denomDF = nVec - t
    ncpVec = nVec * (f * mu)^2/s^2
    critValVec = qf(1 - .05,numDF,denomDF,ncp = 0)
    iRow = 0
    nPoints = 0
    for (n_T in nVec) {
      iRow = iRow + 1
      denomDF = n_T - p
      ncp = n_T * (f * mu)^2/s^2
      critVal = qf(1 - .05,numDF,denomDF,ncp = 0)
      power = 1 - pf(critVal,numDF,denomDF,ncp = ncp)
      powerMat[iRow,iCol] = power
      if ((n_T %% t) == 0) {
        nPoints = nPoints + 1
        plotPointsDF[nPoints,iCol+1] = power
        plotPointsDF[nPoints,1] = n_T
      }
    }
  }
  plotPointsDF = data.frame(plotPointsDF)
  colnames(plotPointsDF) = c("n_T","f_pt025","f_pt05","f_pt1")
  powerMatDF = data.frame(cbind(nVec,powerMat))
  colnames(powerMatDF) = c("n_T", "f_pt025","f_pt05","f_pt1")
  p = ggplot(powerMatDF, aes(x = n_T, y = f_pt025)) +
    geom_line(aes(x = n_T, y = f_pt025), data = powerMatDF, linetype = "solid") + 
    geom_point(aes(x = n_T, y = f_pt025), data = plotPointsDF) + 
    geom_line(aes(x = n_T, y = f_pt05), data = powerMatDF, linetype = "dashed") + 
    geom_point(aes(x = n_T, y = f_pt05), data = plotPointsDF) + 
    geom_line(aes(x = n_T, y = f_pt1), data = powerMatDF, , linetype = "dotted") +
    geom_point(aes(x = n_T, y = f_pt1), data = plotPointsDF) + 
    theme_classic() + ylim(c(0,1)) + xlim(c(0,.95)) + ylab("Lack of fit power") + xlab("Sample size") +
    theme(axis.text.x = element_text(face="bold",size = 10),
          axis.text.y = element_text(face="bold",size = 10), axis.title.x = element_text(face="bold",size = 10),
          axis.title.y = element_text(face="bold",size = 10), title = element_text(face="bold",size = 16)) +
    theme(plot.subtitle = element_text(size = 10)) + theme(plot.title = element_text(size = 14)) +
    scale_x_continuous(breaks = seq(15,90,15)) + 
    scale_y_continuous(breaks = seq(0,1,.1))
  if (iPlot == 1) {
    p = p + labs(title = TeX("(a) $\\sigma = 3$ and $\\tilde{\\mu} = 70$"),
         subtitle = TeX("$f = 0.025$ (solid), $f = 0.05$ (dashed), $f = 0.1$ (dotted)"))
  }
  if (iPlot == 2) {
    p = p + labs(title = TeX("(b) $\\sigma = 3$ and $\\tilde{\\mu} = 52$"),
         subtitle = TeX("$f = 0.025$ (solid), $f = 0.05$ (dashed), $f = 0.1$ (dotted)"))
  }
  if (iPlot == 3) {
    p = p + labs(title = TeX("(c) $\\sigma = 4.5$ and $\\tilde{\\mu} = 70$"),
         subtitle = TeX("$f = 0.025$ (solid), $f = 0.05$ (dashed), $f = 0.1$ (dotted)"))
  }
  if (iPlot == 4) {
    p = p + labs(title = TeX("(d) $\\sigma = 4.5$ and $\\tilde{\\mu} = 52$"),
         subtitle = TeX("$f = 0.025$ (solid), $f = 0.05$ (dashed), $f = 0.1$ (dotted)"))
  }
  
  for (iLine in seq(15,90,15)) {
    p = p + geom_vline(xintercept = iLine, linetype = "solid", linewidth = .1)
  }
  for (iLine in c(0.7,0.8,0.9)) {
    p = p + geom_hline(yintercept = iLine, linetype = "solid", linewidth = .1)
  }
  if (iPlot == 1) {p1 = p}
  if (iPlot == 2) {p2 = p}
  if (iPlot == 3) {p3 = p}
  if (iPlot == 4) {p4 = p}
}
grid.arrange(p1, p2, p3, p4, ncol = 2)
