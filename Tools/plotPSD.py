import matplotlib.pyplot as plt
import numpy as np
import sys

"""
This code is used to plot the psd from a data file (the last time step is used).

Usages:
 (A) python plotPSD.py psd.data
     Reads the radii of all particles in the last time step in the file psd.data and plots the cumulative volume psd.
 (B) python plotPSD.py psd.data psd.csv
     Plots two psd's, the one in the data file and the psd given in the file psd.csv. The csv file has to contain a cumulative number psd, as outputted by the insertion boundary, e.g. "0.1,0\n0.2,0.5\n0.3,1"
 (C) python plotPSD.py psd.data psd.csv orig.csv
     Plots three psd's, the one in the data file and the two psd's given in psd.csv and orig.csv. The idea is that you can compare teh simulation data with the input of the insertion boundary and a original psd (e.g. the psd before cutoff and scaleup)
"""

## read in the data file name
if len(sys.argv) < 2:
    raise ValueError("Specify a data file to read, e.g. python plotPSD.py psd.data")
else:
    fileName = sys.argv[1]
    print("File name: " + fileName)

## open the data file
# read all lines
lines = open(fileName).readlines()
# find line number "no" of the last time step
no = 0
while True:
    # time and number of particles
    t = float(lines[no].split()[1])
    n = int(lines[no].split()[0])
    if no + n + 1 >= len(lines):
        break
    # print("Skipping %r particles at time %.2f" % (n, t))
    no += n + 1
print("Reading %r particles at time %.2f" % (n, t))
# interpret all values as doubles
data = [[float(d) for d in line.split()] for line in lines[no + 1:]]
# read out radii
radius = np.array([dat[6] for dat in data])
# print a few statistics
print("%d particles, d0 %g, dMean %f, d100 %f um"
      % (len(radius), np.min(2*radius), np.mean(2*radius), np.max(2*radius)))

## plot psd of data file
# define figure size in inch
plt.figure(figsize=(12, 6), dpi=80)
# compute number psd
probN, bins = np.histogram(radius, density=True, bins=500)
probN /= sum(probN)
# compute volume psd
r = bins[:-1]
R = bins[1:]
probV = probN * (R**2+r**2) * (R+r)
probV /= sum(probV)
# compute cumulative volume psd
probCV = np.cumsum(probV)
# compute cumulative number psd
probCN = np.cumsum(probN)
# plot cumulative volume psd
plt.plot(2e3*R, probCV,label='simulation')
#plt.plot(r, probCV,label='data')
#plt.plot(R, probCN,label='data CN')
#plt.plot(R, probV/np.max(probV),label='data V')

## read in cumulative number psd from csv
if len(sys.argv) >2:
    csvName = sys.argv[2]
    print("Reading csv file: " + fileName)
    lines = open(csvName).readlines()
    data = [line.strip().split(',') for line in lines]
    radius = np.array([float(line[0]) for line in data])
    # get cumulative number psd
    probCN = np.array([float(line[1]) for line in data])
    # compute number psd
    probN = np.append(probCN[0], np.diff(probCN))
    # compute volume psd
    R = radius
    r = np.append(radius[0],radius[:-1])
    probV = probN * (R**2+r**2) * (R+r)
    probV /= sum(probV)
    # compute cumulative volume psd
    probCV = np.cumsum(probV)
    # plot cumulative volume psd, label as "input"
    plt.plot(2e3*radius,probCV,'g:',label='input')
    #plt.plot(radius,probCN,'r:',label='csv CN')
    #plt.plot(radius,probV/np.max(probV),'k:',label='csv V')
    print("%d particles, d0 %g, dMean %f, d100 %f um (csv)" % (
    len(radius), np.min(2*radius), np.mean(2*radius), np.max(2*radius)))

## read in second cumulative number psd from csv
if len(sys.argv) >3:
    csvName = sys.argv[3]
    print("Reading csv file: " + fileName)
    lines = open(csvName).readlines()
    data = [line.strip().split(',') for line in lines]
    radius = np.array([float(line[0]) for line in data])
    # get cumulative number psd
    probCN = np.array([float(line[1]) for line in data])
    # compute number psd
    probN = np.append(probCN[0], np.diff(probCN))
    # compute volume psd
    R = radius
    r = np.append(radius[0],radius[:-1])
    probV = probN * (R**2+r**2) * (R+r)
    probV /= sum(probV)
    # compute cumulative volume psd
    probCV = np.cumsum(probV)
    # plot cumulative volume psd, label as "original"
    plt.plot(2e3*radius,probCV,'g:',label='original')
    #plt.plot(radius,probCN,'r:',label='csv CN')
    #plt.plot(radius,probV/np.max(probV),'k:',label='csv V')
    print("%d particles, d0 %g, dMean %f, d100 %f um (csv)" % (
        len(radius), np.min(2*radius), np.mean(2*radius), np.max(2*radius)))

## set plot options
plt.ylabel('cumulative volume')
plt.xlabel('d [mm]')
plt.xscale('log')
plt.axis('tight')
plt.legend()

## show and save figure
plt.show()
#plt.savefig('PSD.png')
