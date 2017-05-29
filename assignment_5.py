import numpy as np
import matplotlib.pylab as plt
from scipy.stats import norm
from scipy.stats import chi2
from scipy.interpolate import interp1d

np.random.seed(99)
#domain   = np.linspace(lo,hi,200)


def ComputePoissonAB(n,alpha,beta,lo, hi):
    # Tabulate chisq cdf                                                                                             
    chsqA = []
    chsqB = []
    domain   = np.linspace(lo,hi,26)

    for val in domain:
        chsqA += [chi2.cdf(val, df=(2. * n))]
        chsqB += [chi2.cdf(val, df=(2. * (n + 1)))]
    invChsqA = interp1d(chsqA, domain, bounds_error=False, fill_value=(0.0,51.0))
    invChsqB = interp1d(chsqB, domain, bounds_error=False, fill_value=(0.0,51.0))

    a        = 0.5 * invChsqA(beta)
    b        = 0.5 * invChsqB(1 - alpha)

    if (n == 0 and beta != 0):
        b = - np.log(beta)
        a = 0
        print b

    if (beta == 0.):
        a = 0.

    return a,b



"""

Get the normal sample and test statistic distribution

"""

testStatisticVals  = []

for i in range(0, 10000):
    gaussSample = np.random.normal(loc=100.0, scale=10.0, size=(10))
    xBar        = np.mean(gaussSample)
    lnL0  = len(gaussSample) / 10.0**2  * 100. * (100. - 2. * xBar)
    lnL1  = -len(gaussSample) / 10.0**2  * xBar**2
    testStatisticVals += [lnL0  - lnL1]
    #print testStatisticVals[-1]

hist, bins, _ = plt.hist(testStatisticVals, histtype='step', bins=100, normed=True, color='black')
plt.plot(bins[1:], chi2.pdf(bins[1:],df=1), color='red', lw=2)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel(r"$-2 \ln \Lambda$",fontsize=18)
plt.ylabel(r"PDF",fontsize=18)
plt.xlim(0,10)
plt.ylim(0,2.5)
gaussSample        = np.random.normal(loc=100.0, scale=10.0, size=(10))
gaussSample        = np.random.normal(loc=100.0, scale=10.0, size=(10))
#gaussSample        = np.random.normal(loc=100.0, scale=10.0, size=(10))
#gaussSample        = np.random.normal(loc=100.0, scale=10.0, size=(10))
#gaussSample        = np.random.normal(loc=100.0, scale=10.0, size=(10))


"""

H0: mu = 100, H1: mu != 100

"""

pdfH0 = norm.pdf(gaussSample, loc=100.0, scale=10.0)
xBar  = np.mean(gaussSample)
lnL0  = len(gaussSample) / 10.0**2  * 100. * (100. - 2. * xBar)
lnL1  = -len(gaussSample) / 10.0**2  * xBar**2

print lnL0 - lnL1
dx    = bins[1] - bins[0]
pval  = 0.
i = 0
for x in bins[1:] - dx:
    if (x > lnL0 - lnL1):
        pval += hist[i] * dx
    i += 1

print "xBar is:               " + str(xBar)
print "The pvalue of xbar is: " + str(pval) 
print "Wilk's Thm predicts:   " + str(1 - chi2.cdf(lnL0 - lnL1,df=1))

plt.plot([lnL0  - lnL1, lnL0 - lnL1], [0.,2.5], color='black', linestyle='--', lw=2)
plt.show()


"""

H0: mu = 150, H1: mu != 150

"""

testStatisticVals  = []


for i in range(0, 10000):
    gaussSample = np.random.normal(loc=150.0, scale=10.0, size=(10))
    xBar        = np.mean(gaussSample)
    lnL0  = len(gaussSample) / 10.0**2  * 150. * (150. - 2. * xBar)
    lnL1  = -len(gaussSample) / 10.0**2  * xBar**2
    testStatisticVals += [lnL0  - lnL1]

hist, bins, _ = plt.hist(testStatisticVals, histtype='step', bins=100, normed=True, color='black')

plt.plot(bins[1:], chi2.pdf(bins[1:],df=1), color='red', lw=2)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel(r"$-2 \ln \Lambda$",fontsize=18)
plt.ylabel(r"PDF",fontsize=18)
plt.xlim(0,10)
plt.ylim(0,2.5)

dx    = bins[1] - bins[0]
pval  = 0.
i = 0
for x in bins[1:] - dx:
    if (x > lnL0 - lnL1):
        pval += hist[i] * dx
    i += 1

print "The pvalue of xbar is: " + str(pval)
print "Wilk's Thm predicts:   " + str(1 - chi2.cdf(lnL0 - lnL1,df=1))

plt.plot([lnL0  - lnL1, lnL0 - lnL1], [0.,2.5], color='black', linestyle='--', lw=2)

plt.show()


"""

Compute expected Poisson confidence belt for 10^5 samples.
Keep a running total of how many contain nu = 5.7.

"""

nu         = 5.7
alphaSymm  = 0.05
betaSymm   = 0.05
alphaUp    = 0.1
betaUp     = 0.0
alphaDown  = 0.0
betaDown   = 0.1


hasNuSymm  = 0.
hasNuUp    = 0.
hasNuDown  = 0.

aSymms     = []
bSymms     = []
aUps       = []
bUps       = []
aDowns     = []
bDowns     = []

for i in range(0,100000):
    n            = np.random.poisson(lam=5.7)
    aSymm,bSymm  = ComputePoissonAB(n, alphaSymm, betaSymm, 0., 50.)
    aUp, bUp     = ComputePoissonAB(n, alphaUp, betaUp, 0., 50.)
    aDown, bDown = ComputePoissonAB(n, alphaDown, betaDown, 0., 50.)

    aSymms += [aSymm]
    bSymms += [bSymm]
    aUps   += [aUp]
    bUps   += [bUp]
    aDowns += [aDown]
    bDowns += [bDown]

    if (nu > aSymm and nu < bSymm):
        hasNuSymm += 1
    if (nu < bUp):
        hasNuUp += 1
    if (nu > aDown):
        hasNuDown += 1

    if ((i + 1) % 10 == 0):
        print "Sample " + str(i + 1) + " of 100"

print "The fraction of symmetric intervals containing nu is: "
print float(hasNuSymm)/1e5
print

print "The fraction of upper intervals containing nu is: "
print float(hasNuUp)/1e5
print

print "The fraction of lower intervals containing nu is: "
print float(hasNuDown)/1e5
print

print np.mean(aSymms), np.mean(bSymms)
print np.mean(aUps), np.mean(bUps)
print np.mean(aDowns), np.mean(bDowns)

#print np.arange(0,20)
#print aSymms
plt.hist(aSymms,histtype='step', bins=np.linspace(0,20, 15))
plt.hist(bSymms, histtype='step',bins=np.linspace(0,20,15))
plt.xlabel(r"$\hat{\mu}$", fontsize=18)
plt.ylabel(r"$N$", fontsize=18)
plt.show()

plt.hist(aDowns,histtype='step', bins=np.linspace(0,20, 15))
plt.hist(bUps, histtype='step',bins=np.linspace(0,20,15))
plt.xlabel(r"$\hat{\mu}$", fontsize=18)
plt.ylabel(r"$N$", fontsize=18)
plt.show()
