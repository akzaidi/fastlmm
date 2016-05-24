EstimateLambda = function(fileName)
{
    data = read.delim(fileName, TRUE)
	pvals = data$Pvalue
    LOD2 = median(qchisq(1 - pvals, df = 1))
	return(LOD2 / 0.456)
}

# Tries to create a P-value QQ-plot in -log10(P-value) space and save it to file.
# Gives a warning instead of an error if it fails.
TryMakeQQPlot = function(pvalFileName, outPngFileName)
{
    tryCatch(suppressWarnings(MakeQQPlot(pvalFileName, outPngFileName)),
        error = function(err) warning("can't create QQ plot", call. = FALSE))
}

# Creates a P-value QQ-plot in -log10(P-value) space and saves it to file.
MakeQQPlot = function(pvalFileName, outPngFileName)
{
    pval = read.delim(pvalFileName, TRUE)$Pvalue
    pval[pval > 1] = 1
    
    M = length(pval)
    pNull = (1:M - 0.5)/M
    qNull = -log10(pNull)
    qEmp = -log10(sort(pval))
    
    alphaLevel = 0.05 # significance level for the error bars
    errorBar = CalculateQQErrorBar(M, alphaLevel)
    lowerBar = -log10(errorBar$theoreticalPval - errorBar$betaDown)
    upperBar = -log10(errorBar$theoreticalPval + errorBar$betaUp)
    yBar = -log10(errorBar$theoreticalPval)
    
    png(outPngFileName, width = 800, height = 600, bg = 'white')
    par(mar = c(5, 5.5, 2, 2))
    plot(qNull, qEmp, type = 'p', pch = 19, col = 'blue',
        cex.axis = 1.6, cex.lab = 1.6,
        xlab = expression(-log[10](italic(p))~expected),
        ylab = expression(-log[10](italic(p))~observed),
        xlim = c(0, ceiling(max(qNull))),
        ylim = c(0, ceiling(max(qEmp))))
    lines(yBar, lowerBar, lty = 4, lwd = 2, col = 'dark green')
    lines(yBar, upperBar, lty = 4, lwd = 2, col = 'dark green')
    lines(qNull, qNull, col = "red")
    dev.off()
}

CalculateQQErrorBar = function(M, alphaLevel)
{
    mSeq = 10^(seq(log10(0.5), log10(M - 0.5) + 0.1, 0.1))
    betaQuantiles = data.frame(t(sapply(mSeq,
        function(m) qbeta(c(alphaLevel, 0.5, 1 - alphaLevel), m, M - m))))
    names(betaQuantiles) = c('alpha', 'half', 'oneMinusAlpha')
    
    betaDown = betaQuantiles$half - betaQuantiles$alpha
    betaUp = betaQuantiles$oneMinusAlpha - betaQuantiles$half
    theoreticalPval = mSeq / M
    
    return(list(betaDown = betaDown, betaUp = betaUp, theoreticalPval = theoreticalPval))
}

WritePCCovar = function(covarFileName, pcFileName, T)
{
	covar = read.delim(covarFileName, FALSE)
	pc = read.delim(pcFileName, FALSE)
    if (T == 0)
    {
        pcCovar = covar
    } else
    {
        pcCovar = cbind(covar, pc[, 1:T])
    }
    
    pcCovarFileName = paste('covariates_pc', T, '.txt', sep = '')
	write.table(pcCovar, file = pcCovarFileName, quote = FALSE, sep = '\t',
        row.names = FALSE, col.names = FALSE)
	
    return(pcCovarFileName)
}
