#!/usr/bin/env python3
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

import sys
import os

'''
Short script to combine restart files that are generated in the parallel code, to a single restart file with all
particles and interactions. It is assumed that the restart files are all the same except for the particles and
interactions.

Usage:
- In case of restart files of the form myRestartFile.restart0, myRestartFile.restart1, etc.:
    ./CombineParallelRestartFiles.py myRestartFile.restart
- In case of restart files of the form myRestartFile.restart0.123456, myRestartFile.restart1.123456, etc:
    ./CombineParallelRestartFiles.py myRestartFile.restart 123456
'''


def main():
    # the user gives the base of the restart-file name, we find the restart files of all the cores
    # the argument should include the .restart!
    restart_file_base = sys.argv[1]
    files = sorted([file for file in os.listdir(".") if file.startswith(restart_file_base)])

    # Optionally, the user can give a number s.t. all the files that end with that number will be combined. Necessary
    # if the files are generated with FileType::MULTIPLE_FILES or FileType::MULTIPLE_FILES_PADDED
    if len(sys.argv) > 2:
        restart_file_number = sys.argv[2]
        files = sorted([file for file in files if file.endswith(restart_file_number)])

    print("Combining the following restart files: ")
    print(files)

    # construct the combined file
    out_file = open("combined" + restart_file_base, 'w')

    # The algorithm below is very naive, but works: first construct the header, particles, interactions and footer
    # separately, then write them all to the new file
    header = []
    particles = []
    interactions = []
    footer = []
    footer_written = False

    # Construct the header (everything before the particles) and determine the line number where the particles start
    file0 = open(files[0])
    file0 = file0.readlines()

    index_particles = 0
    while 'Particles' not in file0[index_particles]:
        if 'numberOfProcessors' in file0[index_particles]:
            i = file0[index_particles].find('numberOfProcessors')
            file0[index_particles] = file0[index_particles][:i] + '\n'
        header.append(file0[index_particles])
        index_particles += 1

    # For all files, collect the particle and interaction data. Also collect the footer information once.
    for file in files:
        # Open the file and mold it into a list of lists of all the 'words' in the file
        file0 = open(file)
        file0 = file0.readlines()
        file0_split = [line.split() for line in file0]

        # Collect all particles (we cannot directly write to file, as we need the number of particles)
        index_particles_end = index_particles + int(file0_split[index_particles][1]) + 1
        for i in range(index_particles + 1, index_particles_end):
            particles.append(file0[i])

        # Collect all interactions
        index_interactions_end = index_particles_end + int(file0_split[index_particles_end][1]) + 1
        for i in range(index_particles_end + 1, index_interactions_end):
            interactions.append(file0[i])

        # Collect everything below the interactions, such as h-grid info, chute line
        if not footer_written:
            footer_written = True
            for i in range(index_interactions_end, len(file0_split)):
                footer.append(file0[i])

    # Finally, write everything to the combined restart file
    out_file.write("".join(header))
    out_file.write("Particles " + str(len(particles)) + "\n")
    out_file.write("".join(particles))
    out_file.write("Interactions " + str(len(interactions)) + "\n")
    out_file.write("".join(interactions))
    out_file.write("".join(footer))


if __name__ == '__main__':
    main()
