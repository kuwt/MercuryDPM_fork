#!/usr/bin/env python
#
# Perform Hilbert-Huang transform to the fingering.

import sys
import numpy as np
from scipy import interpolate
from scipy.signal import argrelmax, argrelmin, hilbert


class HHT:
    def __init__(self, Exp):
        self.Exp = Exp
        pars = {}
        with open('%s/Exp.config' % Exp, 'r') as file:
            for line in file:
                info = line.rstrip().split()
                if info:
                    pars[info[0]] = float(info[1])
        self.zrange = pars['width']
        self.radius = pars['radius_large']
        self.final = '%s/Exp.data.%i' % (self.Exp, int(
            pars['timeMax'] / pars['timeStep'] / pars['saveCount']))
        self.front = []
        self.zgrid = np.linspace(-self.zrange, self.zrange,
                                 num=int(self.zrange / self.radius) + 1,
                                 endpoint=True)

        data = []
        with open(self.final, 'r') as file:
            file.readline()
            for line in file:
                info = line.split()
                data.append([float(info[2]), float(info[0])])

        data.sort()
        z = -self.zrange
        grid = []

        for z0, x0 in data:
            if z0 < -self.zrange:
                continue
            elif z - self.radius <= z0 <= z + self.radius:
                grid.append(x0)
            else:
                grid.sort(reverse=True)
                i = 0
                while grid[i] - grid[i + 4] > self.radius * 2:
                    i += 1
                self.front.append(grid[i])
                if z <= self.zrange:
                    grid = [x0]
                    z += self.radius * 2
                else:
                    break
        else:
            grid.sort(reverse=True)
            i = 0
            while grid[i] - grid[i + 4] > self.radius * 2:
                i += 1
            self.front.append(grid[i])

    def stdDev(self, residue1, residue0):
        return sum(i ** 2 for i in (residue0 - residue1)) \
            / sum(i ** 2 for i in residue0)

    def mirrorExtension(self, data):
        return np.concatenate((data[::-1], data, data[::-1]))

    def zgridExtension(self, extrema):
        ext = np.array(self.zgrid)[extrema]
        return np.array(sorted(np.concatenate(
            (-2 * self.zrange - ext, ext, 2 * self.zrange - ext))))

    def sift(self, data, tolerence):
        n = 1
        residue = data.copy()
        while True:
            maxima = argrelmax(residue)[0].tolist()
            minima = argrelmin(residue)[0].tolist()
            if min(len(maxima), len(minima)) <= 1:
                return residue, data - residue
            upper = interpolate.splrep(self.zgridExtension(
                maxima), self.mirrorExtension(residue[maxima]), s=0)
            lower = interpolate.splrep(self.zgridExtension(
                minima), self.mirrorExtension(residue[minima]), s=0)
            mean = (interpolate.splev(self.zgrid, upper) +
                    interpolate.splev(self.zgrid, lower)) / 2
            residue -= mean
            if self.stdDev(residue, residue + mean) < tolerence:
                return residue, data - residue
            n += 1
            if n > 100:
                print "Maximum number of sifting exceeded."
                return residue, data - residue

    def EMD(self, tolerence):
        mode, residue = self.sift(np.array(self.front), tolerence)
        self.residues = [residue.tolist()]
        self.modes = [mode.tolist()]
        while min(len(argrelmax(residue)[0]), len(argrelmin(residue)[0])) > 1:
            if len(self.modes) > 20:
                print "Maximum number of modes exceeded."
                break
            mode, residue = self.sift(residue, tolerence)
            self.residues.append(residue.tolist())
            self.modes.append(mode.tolist())

    def Hilbert(self):
        self.frequencies = []
        for mode in self.modes:
            analytic_signal = hilbert(mode)
            phase = np.unwrap(np.angle(analytic_signal))
            frequency = np.diff(phase) / (2 * np.pi * self.radius)
            self.frequencies.append([0] + frequency.tolist())

    def write(self):
        with open(self.Exp + '/front.data', 'w') as file:
            for info in zip(self.zgrid, self.front):
                file.write(' '.join(map(str, info)) + '\n')
        with open(self.Exp + '/EMD.data', 'w') as file:
            for info in zip(self.zgrid, *self.modes):
                file.write(' '.join(map(str, info)) + '\n')
        with open(self.Exp + '/Hilbert.data', 'w') as file:
            for info in zip(self.zgrid, *self.frequencies):
                file.write(' '.join(map(str, info)) + '\n')


def main(Exp):
    data = HHT(Exp)
    data.EMD(1e-5)
    data.Hilbert()
    data.write()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: %s Exp\n' % sys.argv[0]
        exit(1)
    main(sys.argv[1])
