#!/usr/bin/env python3
#
# Perform individual large particles tracking.

import sys
import math
import random


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
        self.tracking = []
        self.Lagrangian = [[] for _ in xrange(10)]

    def retrieveMax(self, data):
        data.sort(reverse=True)
        i = 0
        while data[i] - data[i + 10] > 2 * self.radius:
            i += 1
        return data[i]

    def readLast(self):
        particles = []
        data = []
        with open('%s/front.data' % self.Exp) as file:
            for line in file:
                data.append(float(line.split()[1]))
        front = max(data)
        with open('%s/Exp.data.%i' % (self.Exp, self.last)) as file:
            file.readline()
            for line in file:
                info = line.split()
                if info[-1] == '1' and float(info[0]) > front * 0.99:
                    particles.append(map(float, info[:3]) + [info[6]])
        return particles

    def read(self, i):
        output = []
        data = []
        with open('%s/Exp.data.%i' % (self.Exp, i)) as file:
            file.readline()
            for line in file:
                info = line.split()
                data.append(float(info[0]))
                if info[6] in self.tracking:
                    x, y, z = map(float, info[:3])
                    output.append([x, y + x * self.angle, z])
        front = self.retrieveMax(data)
        for i, (x, y, z) in enumerate(output):
            self.Lagrangian[i].append([x - front, y, z])
        return output

    def write(self, i, particles):
        with open('%s/track.data.%i' % (self.Exp, i), 'w') as file:
            for x, y, z in particles:
                file.write('%f %f %f\n' % (x, y, z))

    def track(self):
        particles = self.readLast()
        output = []
        for particle in random.sample(particles, 10):
            self.tracking.append(particle[3])
        print 'Particles selected.'
        for i in xrange(self.last + 1):
            output = self.read(i)
            self.write(i, output)
            print '%i completed.' % i
        for i in xrange(10):
            with open('%s/Lagrangian.data.%i' % (self.Exp, i), 'w') as file:
                for x, y, z in self.Lagrangian[i]:
                    file.write('%f %f %f\n' % (x, y, z))


def main(Exp):
    flow = Profile(Exp)
    flow.track()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: %s Exp\n' % sys.argv[0]
        exit(1)
    main(sys.argv[1])
