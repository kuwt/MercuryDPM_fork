#!/usr/bin/env python
#Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
#For the list of developers, see <http://www.MercuryDPM.org/Team>.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name MercuryDPM nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os.path
import sys

print("Merging parallel generated restart files.")

simulation_name = sys.argv[1:][0]
print simulation_name
print ("Simulation name is " + simulation_name)

#Open the first restart file. This is always name.restart0
root_filename = simulation_name + '.restart0'
if not os.path.isfile(root_filename):
    print("Error: Cannot find " + root_filename)
    sys.exit(0)
f = open(root_filename, 'r+')

#Find the number of cores
found = False
while not found:
    line = f.readline()
    tokens = line.split()
    for i, j in enumerate(tokens):
        if j == "numberOfProcessors":
            nProc = int(tokens[i+1])
            print ("Simulation has run on " + str(nProc) + " number of cores")
            found = True
        if j == "":
            print("Error: Can not find the number of cores")
f.close()

#Check if other files are present
for i in range(1,nProc):
    filename = simulation_name + ".restart" + str(i)
    if not os.path.isfile(filename):
        print("Error: Cannot find " + filename)
        sys.exit(0)
print("File check: All restart files have been found")

#Find the number of particles and number of interactions in each file
nParticles = []
nInteractions = []
for i in range(0,nProc):
    filename = simulation_name + ".restart" + str(i)
    f_temp = open(filename, 'r')
    found_particles = False
    found_interactions = False
    while not (found_particles and found_interactions):
        line = f_temp.readline()
        tokens = line.split()
        if len(tokens) == 2:
            if tokens[0] == 'Particles':
                nParticles.append(int(tokens[1]))
                found_particles = True
            if tokens[0] == 'Interactions':
                nInteractions.append(int(tokens[1]))
                found_interactions = True
    f_temp.close()
nParticlesTotal = sum(nParticles)
nInteractionsTotal = sum(nInteractions)
print(str(nParticles))
print(str(nParticlesTotal))
print("Total number of particles: " + str(nParticlesTotal))
print("Total number of interactions: "+ str(nInteractionsTotal))

#Create new file to write
filename_new = simulation_name + ".restart"
print("Creating: " + filename_new)
open(filename_new, 'w').close() #Empty the file if it exists
f_new = open(filename_new,'a')
f = open(root_filename, 'r+')

#Write generic data from the root restart file to the new file
finished = False
f.seek(0)
while not finished:
    line = f.readline()
    tokens = line.split()
    if not (tokens[0] == 'Particles'):
        f_new.write(line)
    else: 
        finished = True

#Write particle data from all files into the new file
f_new.write("Particles " + str(nParticlesTotal) + "\n")
for i in range (0,nProc):
    filename = simulation_name + ".restart" + str(i)
    f_temp = open(filename, 'r')
    found_particles = False
    found_interactions = False
    while not (found_particles and found_interactions):
        line = f_temp.readline()
        tokens = line.split()

        #Check if the current line has reached the interactions
        if (tokens[0] == 'Interactions'):
            found_interactions = True

        #Copy the line from the temp file to the new file
        if (found_particles == True and found_interactions == False):
            f_new.write(line)

        #Check if we have reached the particle lines
        if tokens[0] == 'Particles':
            found_particles = True
    f_temp.close()

#Write interaction data from all files into the new file
f_new.write("Interactions " + str(nInteractionsTotal) + "\n")
for i in range (0,nProc):
    filename = simulation_name + ".restart" + str(i)
    f_temp = open(filename, 'r')
    found_end = False 
    found_interactions = False
    while not (found_particles and found_end):
        line = f_temp.readline()
        tokens = line.split()
        if (len(tokens) == 0):
            found_end = True
        else:
            #Check if the current line has reached the interactions
            if (tokens[0] == 'NUM_BUCKETS' or tokens[0] == ""):
                found_end = True

            #Copy the line from the temp file to the new file
            if (found_interactions == True and found_end == False):
                f_new.write(line)

            #Check if we have reached the particle lines
            if tokens[0] == 'Interactions':
                found_interactions = True
    f_temp.close()

#Finally write Hgrid
finished = False
f.seek(0)
while not finished:
    line = f.readline()
    tokens = line.split()
    if (tokens[0] == 'NUM_BUCKETS'):
        f_new.write(line)
        finished = True

#Close all files
f.close()
f_new.close()
print("Finished writing file: " + filename_new)






























