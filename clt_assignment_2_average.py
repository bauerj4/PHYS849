import random
import thread
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np

random.seed(1337)

"""

Get k random distinct elements and 
return an array of length n, the joint
k-sum distribution

"""

def gaussPDF(x,mu,sigma):
    val = np.exp(-(x-mu)**2 / (2. * sigma**2))
    val /= (2. * np.pi)**0.5 * sigma
    return val

vGaussPDF = np.vectorize(gaussPDF, excluded=['mu','sigma'])

def selfSampleSum(n,k,arr):
    sumSample = []

    for i in range(n):
        theseIndices = dict([])
        thisSum = 0
        i = 0
        while (i < k):
            idx = random.randint(0, len(arr) - 1)
            if (not theseIndices.has_key(idx)):
                thisSum += arr[idx]
                theseIndices[idx] = 1
                i += 1
        sumSample += [thisSum]
    return sumSample

unifArr = [random.uniform(0,1) for i in range(10000)]

"""

The uniform distribution histogram

"""

bins = np.linspace(0,1,25)
dx   = bins[1] - bins[0]

meanUnif  = np.mean(unifArr)
stdUnif   = np.std( unifArr)

print meanUnif, stdUnif

plt.hist(unifArr,histtype='step',bins=bins)
plt.plot([0,1],[10000./24., 10000./24.], color='black', lw=2)
plt.xlabel("$x$",fontsize=20)
plt.ylabel("$N$ (raw counts)", fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

ax = plt.gca()
ax.annotate(r'$\mu$ = {:f}, $\sigma$ = {:f}'.format(meanUnif, stdUnif), xy=(0.1,360),fontsize=18)

plt.show()

"""

Sums of elements

"""
stdGaussSpace          = np.linspace(0,1,1000)

twoWaySumElements      = selfSampleSum(10000,2,   unifArr)
tenWaySumElements      = selfSampleSum(10000,10,  unifArr)
hundredWaySumElements  = selfSampleSum(10000,100, unifArr)
thousandWaySumElements = selfSampleSum(10000,1000,unifArr)


twoWayAverage          = np.array(twoWaySumElements)/2.
tenWayAverage          = np.array(tenWaySumElements)/10.
hundredWayAverage      = np.array(hundredWaySumElements)/100.
thousandWayAverage     = np.array(thousandWaySumElements)/1000.

plt.figure(figsize=(16,6))

twoAverageMean       = np.mean(twoWayAverage)
tenAverageMean       = np.mean(tenWayAverage)
hundredAverageMean   = np.mean(hundredWayAverage)
thousandAverageMean  = np.mean(thousandWayAverage)

twoAverageStd        = np.std(twoWayAverage) #/ np.sqrt(2.)
tenAverageStd        = np.std(tenWayAverage) #/ np.sqrt(10.)
hundredAverageStd    = np.std(hundredWayAverage)# / np.sqrt(100.)
thousandAverageStd   = np.std(thousandWayAverage)# / np.sqrt(1000.)

print twoAverageMean, twoAverageStd
print tenAverageMean, tenAverageStd
print hundredAverageMean, hundredAverageStd
print thousandAverageMean, thousandAverageStd

a2     = 1./(np.sqrt(2. * np.pi) * twoAverageStd * dx)
a10    = 1./(np.sqrt(2. * np.pi) * tenAverageStd * dx)
a100   = 1./(np.sqrt(2. * np.pi) * hundredAverageStd * dx)
a1000  = 1./(np.sqrt(2. * np.pi) * thousandAverageStd * dx)

print a2, a10, a100, a1000

twoAverageGauss       = vGaussPDF(stdGaussSpace, twoAverageMean, twoAverageStd)
tenAverageGauss       = vGaussPDF(stdGaussSpace, tenAverageMean, tenAverageStd)
hundredAverageGauss   = vGaussPDF(stdGaussSpace, hundredAverageMean, hundredAverageStd)
thousandAverageGauss  = vGaussPDF(stdGaussSpace, thousandAverageMean, thousandAverageStd)

gs = gridspec.GridSpec(1,2)

ax1 = plt.subplot(gs[0])
hist1, _, _ = plt.hist(twoWayAverage,      histtype='step', color='blue',  normed=True)
hist2, _, _ = plt.hist(tenWayAverage,      histtype='step', color='green', normed=True)
hist3, _, _ = plt.hist(hundredWayAverage,  histtype='step', color='red',   normed=True)
hist4, _, _ = plt.hist(thousandWayAverage, histtype='step', color='teal',  normed=True)

plt.plot(stdGaussSpace, twoAverageGauss, color='blue', lw=2,ls='--')
plt.plot(stdGaussSpace, tenAverageGauss, color='green', lw=2,ls='--')
plt.plot(stdGaussSpace, hundredAverageGauss, color='red', lw=2,ls='--')
plt.plot(stdGaussSpace, thousandAverageGauss, color='teal', lw=2,ls='--')

plt.xlabel(r"$\bar{x}$",fontsize=20)
plt.ylabel("$N$ (normalized counts)", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


ax2 = plt.subplot(gs[1])
ax2.set_yscale('log')

plt.ylim(1e-3,75)

hist1, _, _ = plt.hist(twoWayAverage,      histtype='step', color='blue',  normed=True)
hist2, _, _ = plt.hist(tenWayAverage,      histtype='step', color='green', normed=True)
hist3, _, _ = plt.hist(hundredWayAverage,  histtype='step', color='red',   normed=True)
hist4, _, _ = plt.hist(thousandWayAverage, histtype='step', color='teal',  normed=True)

plt.plot(stdGaussSpace, twoAverageGauss, color='blue', lw=2,ls='--')
plt.plot(stdGaussSpace, tenAverageGauss, color='green', lw=2,ls='--')
plt.plot(stdGaussSpace, hundredAverageGauss, color='red', lw=2,ls='--')
plt.plot(stdGaussSpace, thousandAverageGauss, color='teal', lw=2,ls='--')

plt.xlabel(r"$\bar{x}$",fontsize=20)
plt.ylabel("$N$ (normalized counts)", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


plt.show()


plt.figure(figsize=(18,10))
baseStd  = (1./12.)**0.5
baseMu   = 0.5
gs = gridspec.GridSpec(2,2)



ax1 = plt.subplot(gs[0])
thisMu     = 2. * baseMu
thisStd    = 2. * baseStd / np.sqrt(2)  # Var = N thisVar
gaussRange = np.linspace(thisMu - 3. * thisStd, thisMu + 3. * thisStd, 101)
gauss      = vGaussPDF(gaussRange,thisMu,thisStd)
plt.hist(twoWaySumElements, histtype='step',normed='true',bins=25)
plt.plot(gaussRange, gauss, color='black', lw=2)
plt.xlabel("$x$",fontsize=20)
plt.ylabel("$N$ (normalized counts)", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.xlim(gaussRange[0], gaussRange[-1])
plt.ylim(0, 1.25 * np.amax(gauss))
ax1.annotate(r'$N_p$ = 2, $\mu$ = {:f}, $\sigma$ = {:f}'.format(np.mean(twoWaySumElements), \
                                                                 np.std(twoWaySumElements)),\
             xy=(thisMu-2.5 *thisStd,1.1*np.amax(gauss)),fontsize=18)


ax2 = plt.subplot(gs[1])
thisMu     = 10. * baseMu
thisStd    = 10. * baseStd / np.sqrt(10.)
gaussRange = np.linspace(thisMu - 3. * thisStd, thisMu + 3. * thisStd, 101)
gauss      = vGaussPDF(gaussRange,thisMu,thisStd)
plt.hist(tenWaySumElements, histtype='step',normed='true',bins=25)
plt.plot(gaussRange,gauss, color='black', lw=2)
plt.xlabel("$x$",fontsize=20)
plt.ylabel("$N$ (normalized counts)", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

ax2.annotate(r'$N_p$ = 10, $\mu$ = {:f}, $\sigma$ = {:f}'.format(np.mean(tenWaySumElements), \
                                                                 np.std(tenWaySumElements)),\
             xy=(thisMu-2.5 *thisStd,1.1*np.amax(gauss)),fontsize=18)

plt.xlim(gaussRange[0], gaussRange[-1])
plt.ylim(0, 1.25 * np.amax(gauss))


ax3 = plt.subplot(gs[2])
thisMu     = 100. * baseMu
thisStd    = 100. * baseStd / np.sqrt(100.)
gaussRange = np.linspace(thisMu - 3. * thisStd, thisMu + 3. * thisStd, 101)
gauss      = vGaussPDF(gaussRange,thisMu,thisStd)
plt.hist(hundredWaySumElements, histtype='step',normed='true',bins=25)
plt.plot(gaussRange,gauss, color='black', lw=2)
plt.xlabel("$x$",fontsize=20)
plt.ylabel("$N$ (normalized counts)", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

ax3.annotate(r'$N_p$ = 100, $\mu$ = {:f}, $\sigma$ = {:f}'.format(np.mean(hundredWaySumElements), \
                                                                 np.std(hundredWaySumElements)),\
             xy=(thisMu-2.5 *thisStd,1.1*np.amax(gauss)),fontsize=18)

plt.xlim(gaussRange[0], gaussRange[-1])
plt.ylim(0, 1.25 * np.amax(gauss))


ax4 = plt.subplot(gs[3])
thisMu     = 1000. * baseMu
thisStd    = 1000. * baseStd / np.sqrt(1000.)
gaussRange = np.linspace(thisMu - 3. * thisStd, thisMu + 3. * thisStd, 101)
gauss      = vGaussPDF(gaussRange,thisMu,thisStd)
plt.hist(thousandWaySumElements, histtype='step',normed='true',bins=25)
plt.plot(gaussRange, gauss, color='black', lw=2)
plt.xlabel("$x$",fontsize=20)
plt.ylabel("$N$ (normalized counts)", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

ax4.annotate(r'$N_p$ = 1000, $\mu$ = {:f}, $\sigma$ = {:f}'.format(np.mean(thousandWaySumElements), \
                                                                 np.std(thousandWaySumElements)),\
             xy=(thisMu-2.5 *thisStd,1.1*np.amax(gauss)),fontsize=18)

plt.xlim(gaussRange[0], gaussRange[-1])
plt.ylim(0, 1.25 * np.amax(gauss))

plt.show()
#plt.hist(unifArr, histtype='step')
#plt.show()
