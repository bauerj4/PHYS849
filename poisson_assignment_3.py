from scipy.stats import poisson
from scipy.stats import binom

mu  = 100.0
std = 10.0

poissonThreeStd  = 0.
poissonFiveStd   = 0.
anyBinGTThreeStd = 0.
anyBinGTFiveStd  = 0.

for i in range(130,200):
    poissonThreeStd += poisson.pmf(i,mu)

for i in range(150, 200):
    poissonFiveStd += poisson.pmf(i,mu)

for i in range(0, 21):
    anyBinGTThreeStd += 1. - binom.cdf(i,20,poissonThreeStd)

for i in range(0, 21):
    anyBinGTFiveStd += 1. - binom.cdf(i,20,poissonFiveStd)


print "P[bin count > 130] = " + str(poissonThreeStd)
print "P[bin count > 150] = " + str(poissonFiveStd)
print "P[any bin > 130]   = " + str(anyBinGTThreeStd)
print "P[any bin > 150]   = " + str(anyBinGTFiveStd)

