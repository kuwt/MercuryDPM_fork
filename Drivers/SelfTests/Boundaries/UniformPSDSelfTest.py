import matplotlib.pyplot as plt
import numpy as np

# open data file
fileName = "UniformPSDSelfTest.data"
print("Reading "+fileName)
lines = open(fileName).readlines()
data = [[float(d) for d in line.split()] for line in lines[1:]]
# read out size and species
radius = np.array([dat[6] for dat in data])
print("%d particles, radius: min %g, mean %f, max %f" % (len(radius), np.min(radius), np.mean(radius), np.max(radius)))

# plot psd of species 0
probN, bins = np.histogram(radius, density=True, bins=30)
r = bins[:-1]
R = bins[1:]
probV = probN * (R**2+r**2) * (R+r)
probV /= sum(probV)
plt.plot(bins[:-1], probV,'-',label='data')
plt.plot(bins[:-1], bins[:-1]**3/np.sum(bins[:-1]**3)/np.sum(probV),'-',label='analytic')
# set plot options
plt.ylabel('vpsd')
plt.xlabel('radius')
plt.legend()
plt.show()
plt.savefig('UniformPSDSelfTest.png',dpi=300)
