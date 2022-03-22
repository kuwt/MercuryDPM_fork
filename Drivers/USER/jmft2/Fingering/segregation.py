#!/usr/bin/env python3
#
# Generate velocity profile of large particles at different height.

import sys
import math
from collections import defaultdict


class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __add__(self, vector):
        return Vector(self.x + vector.x, self.y + vector.y)

    def __div__(self, scaler):
        return Vector(self.x / scaler, self.y / scaler)

    def rot(self, angle):
        self.x = self.x * math.cos(angle) - self.y * math.sin(angle)
        self.y = self.x * math.sin(angle) + self.y * math.cos(angle)
        return self

    def rescale(self):
        length = math.sqrt(self.x**2 + self.y**2)
        if length:
            new = self / (100 * length)
            return new.x, new.y, length
        else:
            return 0, 0, 0


class Particle:
    def __init__(self, info, slope):
        self.x = float(info[0])
        self.y = float(info[1]) + self.x * slope
        self.z = float(info[2])
        self.vx = float(info[3])
        self.vy = float(info[4])
        self.vz = float(info[5])
        self.species = int(info[-1])


class Grid:
    def __init__(self, particles, size):
        self.size = size
        self.grid = defaultdict(list)
        for particle in particles:
            self.addParticle(particle)

    def addParticle(self, particle):
        x = int(round(particle.x / self.size))
        y = int(round(particle.y / self.size))
        self.grid[(x, y)].append(particle)

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
        self.size = self.radius
        self.angle = pars['alpha']
        self.x = int(pars['length'] * 20 / self.size)
        self.y = int(pars['height'] / self.size)

    def read(self, i):
        particles = []
        with open('%s/Exp.data.%i' % (self.Exp, i)) as file:
            file.readline()
            for line in file:
                info = line.split()
                particles.append(Particle(info, math.tan(self.angle)))
        return particles

    def flowProfile(self, grid):
        profile = defaultdict(lambda: ['NaN'])
        for (x, y), particles in grid:
            if len(particles) < 5:
                continue
            N0 = N1 = 0
            v0 = v1 = Vector(0, 0)
            for particle in particles:
                if particle.species == 0:
                    N0 += 1
                    v0 += Vector(particle.vx, particle.vy)
                else:
                    N1 += 1
                    v1 += Vector(particle.vx, particle.vy)
            ratio = float(N1) / (N0 + N1)
            if N0 < 5:
                v0 = Vector(0, 0)
            else:
                v0 /= N0
                v0.rot(self.angle)
            if N1 < 5:
                v1 = Vector(0, 0)
            else:
                v1 /= N1
                v1.rot(self.angle)
            profile[(x, y)] = ratio, (v0.rescale(), v1.rescale())
        return profile

    def write(self, i, profile):
        with open('%s/segregation.data.%i' % (self.Exp, i), 'w') as file:
            for (x, y), (_, (v0, v1)) in profile.iteritems():
                file.write('%f %f ' % (x * self.size, y * self.size))
                file.write(' '.join(map(str, v0 + v1)) + '\n')
        with open('%s/ratio.data.%i' % (self.Exp, i), 'w') as file:
            for y in xrange(self.y + 1):
                for x in xrange(self.x + 1):
                    ratio = profile[(x, y)][0]
                    file.write(str(ratio) + ' ')
                file.write('\n')

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
