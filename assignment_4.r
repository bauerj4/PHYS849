library(ggplot2)
library(graphics)

set.seed(1338)
sampleSize = 1000
nSamples   = 1000


# Binary search for optimal tau

#ChiSqOptimize <- function(bins,counts, tauGuess, domainWidth, bmax, niter){
#    lo        <- tauGuess - domainWidth/2.
#    hi        <- tauGuess + domainWidth/2.
#
#    thisChisq <- 0.
#    for (i in 1:niter){
#        hiChisq   <- DecayChisq(bins,counts,hi,bmax)
#        loChisq   <- DecayChisq(bins,counts,lo,bmax)
#        thisChisq <- DecayChisq(bins,counts,(hi - lo)/2.,bmax)
#
#        if (thisChisq < hiChisq & thisChisq > loChisq){
#            hi <- (hi - lo)/2.
#        }
#        
#    }
#    
#}



# Decay ChiSq

DecayChiSq <- function(tau,bins,counts,bmax){
    dx        <- bins[2] - bins[1]
    residuals <- (counts - sampleSize * exp(-bins/tau) * dx/tau)
                  #(1./tau * (1. - exp(-bmax)) * exp(-bins/tau)) * dx)
    chisq     <- residuals^2 / (sampleSize * exp(-bins/tau) * dx/tau)
                                        #(nSamples * (1./tau * (1. - exp(-bmax)) * exp(-bins/tau)) * dx)
                                        #sd(residuals)^2
    #cat(sum(residuals), sum(counts))
    return((sum(chisq)))
}

# Evaluate the MLE for tau

DecayEstimator <- function(sample){
    return(mean(sample))
}


# Draw a sample from a normal distribution ensuring
# positive values

GaussSample <- function(bmax,n){
    sample <- c()
    for (i in 1:n){
        normsample = rnorm(1, mean=1, sd=0.1)
        if (normsample < 0 || normsample > bmax){
            i = i - 1
        } else {
            sample <- c(sample, normsample)
        }
    }
    return(sample)
}

# Draw a sample from an exponential distribution
# via inversion sampling.  Use truncated normalization.

ExpSample <- function(tau,bmax,n){
    sample <- c()
    norm   <- 1./tau * (1. - exp(-bmax))
    for (i in 1:as.integer(n)){
        unifsample <- runif(1)
        #expsample  <- -tau * log(1./norm * unifsample)
        expsample <- -log(1 - unifsample)*tau
        if (expsample > bmax || expsample < 0){
            i = i - 1
        } else {
            sample     <- c(sample,expsample)
        }
    }
    #print(mean(sample))
    return(sample)
}


expMaxRange   <- 6.5 # Exponential truncation
expDecay      <- 1.0 # Exponential decay constant
expBreaks     <- seq(0.,expMaxRange,length.out=20)
expSamples    <- matrix(,nrow=nSamples,ncol=(length(expBreaks)-1)) # Store the results
expBins       <- c()
decayMLE      <- c()
chisqExp      <- c()
optimTauExp   <- c()
optimChisqExp <- c()

normMaxRange <- 2.
normBreaks   <- seq(0., normMaxRange,length.out=20)
normSamples  <- matrix(,nrow=nSamples,ncol=(length(expBreaks) - 1))
normBins     <- c()
decayNormMLE <- c()
chisqNorm    <- c()


for (i in 1:nSamples){
    exponential   <- ExpSample(expDecay,expMaxRange, sampleSize)
    expNorm       <- 1./sampleSize * 1./expDecay #* (1 - exp(-expMaxRange))
    expHist       <- hist(exponential, breaks=expBreaks)
    expCounts     <- expHist$counts * expNorm / (expHist$mids[2] - expHist$mids[1])
    expBins       <- expHist$mids
    decayMLE      <- c(decayMLE, DecayEstimator(exponential))
    chisqExp      <- c(chisqExp, DecayChiSq(expDecay,expBins, expHist$counts, expMaxRange))

    #print(optim(1.5,DecayChiSq,upper=1.25, lower=0.75, method="SANN", expBins, expHist$counts, expMaxRange))
    
    optimExp      <- optim(1.0,DecayChiSq,method="Nelder-Mead",gr=NULL, expBins, expHist$counts, expMaxRange)
    optimTauExp   <- c(optimTauExp, optimExp$par)
    optimChisqExp <- c(optimChisqExp, optimExp$value)
    
    gauss         <- GaussSample(normMaxRange, sampleSize)
    gaussNorm     <- 1./(0.1 * (3.14159 * 2.)^0.5)
    normHist      <- hist(gauss,breaks=normBreaks)
    normBins      <- normHist$mids
    normCounts    <- normHist$counts / (sampleSize * (normHist$mids[2] - normHist$mids[1]))
    decayNormMLE  <- c(decayNormMLE, DecayEstimator(gauss))
    chisqNorm     <- c(chisqNorm, DecayChiSq(expDecay,normBins, normHist$counts, expMaxRange))
    
    for (j in 1:length(expBins)){
        expSamples[i,j] <- expCounts[j]
    }

    for (j in 1:length(normBins)){
        normSamples[i,j] <- normCounts[j]
    }
    
    if (i %% 10 == 0 || i == nSamples ){
        cat("Drew sample ", i,'\n')
    }
}


#print(chisqExp)
print(chisqNorm)
print(mean(chisqExp))
expMeanCounts  <- c()
expStdCounts   <- c()
normMeanCounts <- c()
normStdCounts  <- c()
for (i in 1:length(expBins)){
    expMeanCounts <- c(expMeanCounts, mean(expSamples[,i]))
    expStdCounts  <- c(expStdCounts, sd(expSamples[,i]))
}

for (i in 1:length(normBins)){
    normMeanCounts <- c(normMeanCounts, mean(normSamples[,i]))
    normStdCounts  <- c(normStdCounts, sd(normSamples[,i]))
}


print(sum(expMeanCounts))

dx               <- expBins[2] - expBins[1]
decayMLEDens     <- density(decayMLE)
decayMLEFn       <- approxfun(decayMLEDens$x, decayMLEDens$y, rule=2)
normDecayMLEDens <- density(decayNormMLE)
normDecayMLEFn   <- density(normDecayMLEDens$x, normDecayMLEDens$y, rule=2)
#chisqExp     <- DecayChiSq(decayMLE,expDecay)
#print(chisqExp)
chisqExpDens  <- density(chisqExp, bw=0.5)
chisqNormDens <- density(chisqNorm)
chisq18       <- dchisq(seq(0,150,length.out=150), 18)


optimTauExpDens   <- density(optimTauExp, bw=0.5)
optimChisqExpDens <- density(optimChisqExp, bw=0.5)
optimTauExpFn     <- approxfun(optimTauExpDens$x, optimTauExpDens$y, rule=2)
optimChisqExpFn   <- approxfun(optimChisqExpDens$x, optimChisqExpDens$y, rule=2)

pvalueFirst  <- integrate(decayMLEFn, decayMLE[1], expMaxRange)

cat(decayMLE[1], "\n")
print(pvalueFirst)

setEPS()
postscript("expsampleR.eps")
plot(expBins,expMeanCounts)
lines(expBins,   exp(-expBins/expDecay))
lines(expBins, expMeanCounts + 2.*expStdCounts,col='red')
lines(expBins, expMeanCounts - 2.*expStdCounts,col='red')


postscript("expoptimizedtau.eps")
plot(optimTauExpDens)

postscript("expoptimizedchisq.eps")
plot(optimChisqExpDens)


postscript("tauestimatorhist.eps")
#hist(decayMLE)
plot(density(decayMLE),main="",xlab="",ylab="")
lines(optimTauExpDens$x, optimTauExpDens$y, lty=2)
title(main="",xlab=expression(tau),ylab=expression(f(tau)), cex.main=1.4, cex.lab=1.4)


postscript("expchisq.eps")
#hist(chisqExp)
plot(chisqExpDens,ylim=c(0.,0.09))
#lines(chisqNormDens)
lines(optimChisqExpDens$x, optimChisqExpDens$y, lty=2)
lines(seq(0,150,length.out=150),chisq18, col='red')


postscript("normtauestimatorhist.eps")
plot(normDecayMLEFn,main="",xlab="",ylab="")
title(main="",xlab=expression(tau),ylab=expression(f(tau)), cex.main=1.4, cex.lab=1.4)


postscript("normsampleR.eps")
plot(normBins,normMeanCounts)
lines(normBins,   exp(-normBins/expDecay))
lines(normBins, normMeanCounts + 2.*normStdCounts,col='red')
lines(normBins, normMeanCounts - 2.*normStdCounts,col='red')


postscript("normchisq.eps")
plot(chisqNormDens)
