#!/usr/bin/env python3
#
# Generate height and velocity profile of the flow.
#
# TODO(bj268): Integrate the code to kernel (DPMBase.cc)
# since most of the time is spent on reading data files.

import sys
import math
from collections import defaultdict


class Particle:
    def __init__(self, x, y, z, vx, vy, vz):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz


class Grid:
    def __init__(self, particles, size):
        self.size = size
        self.grid = defaultdict(list)
        for particle in particles:
            self.addParticle(particle)

    def addParticle(self, particle):
        x = int(round(particle.x / self.size))
        z = int(round(particle.z / self.size))
        self.grid[(x, z)].append(particle)

    def __iter__(self):
        return self.grid.iteritems()


class Profile:
    def __init__(self, Exp):
        pars = {}
        self.Exp = Exp
        with open('%s/Exp.config' % self.Exp, 'r') as file:
            for line in file:
                info = line.rstrip().split()
                if info:
                    pars[info[0]] = float(info[1])
        self.last = int(pars['timeMax'] / pars['timeStep'] / pars['saveCount'])
        self.radius = pars['radius_large']
        self.size = self.radius * 4
        self.angle = math.tan(pars['alpha'])
        self.x = int(pars['length'] * 20 / self.size)
        self.z = int(pars['width'] * 1.2 / self.size)

    def read(self, i):
        particles = []
        with open('%s/Exp.data.%i' % (self.Exp, i)) as file:
            file.readline()
            for line in file:
                info = line.split()
                particles.append(Particle(*map(float, info[:6])))
        return particles

    def flowProfile(self, grid):
        # TODO(bj268): Use depth-averaged velocity instead and
        # improve height determination if required.
        profile = defaultdict(lambda: ('NaN', 'NaN'))
        for (x, z), particles in grid:
            if len(particles) < 5:
                continue
            cell = []
            for particle in particles:
                cell.append((particle.y + particle.x * self.angle,
                             particle.vx**2 + particle.vy**2 + particle.vz**2))
            cell.sort()
            middle = cell[len(cell) / 2 - 2:len(cell) / 2 + 3]
            profile[(x, z)] = (cell[-1][0] + self.radius,
                               math.sqrt(sum(v2 for h, v2 in middle) / 5))
        return profile

    def write(self, i, profile):
        with open('%s/height.data.%i' % (self.Exp, i), 'w') as hfile, \
                open('%s/speed.data.%i' % (self.Exp, i), 'w') as vfile:
            for z in xrange(-self.z, self.z + 1):
                for x in xrange(self.x + 1):
                    h, v = profile[(x, z)]
                    hfile.write(str(h) + ' ')
                    vfile.write(str(v) + ' ')
                hfile.write('\n')
                vfile.write('\n')

    def analyze(self):
        for i in xrange(self.last + 1):
            particles = self.read(i)
            grid = Grid(particles, self.size)
            self.write(i, self.flowProfile(grid))
            print '%i completed.' % i


def main(Exp):
    flow = Profile(Exp)
    flow.analyze()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: %s Exp\n' % sys.argv[0]
        exit(1)
    main(sys.argv[1])
