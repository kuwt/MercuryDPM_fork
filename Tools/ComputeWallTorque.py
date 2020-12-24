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
from math import sqrt

'''
Short script to compute torque on walls from restart files.

Usage:
- To compute torque on wall from the file name.restart:
    ./ConvertRestartToData.py name.restart wall_id
- To compute torque on a range of walls from the file name.restart:
    ./ConvertRestartToData.py name.restart wall_id_begin wall_id_end
'''


def main():

    # the user gives the base of the restart-file name
    if (len(sys.argv)<3):
        raise(Exception("This code needs at least two arguments, e.g. \n\t./ConvertRestartToData.py name.restart wall_id [wall_id_end]"))

    wall_id_begin = int(sys.argv[2])
    if (len(sys.argv)>3):
        wall_id_end = int(sys.argv[3])
    else:
        wall_id_end = wall_id_begin

    # name data and fstat
    restart_file_name = sys.argv[1]
    data_file_name = restart_file_name.replace(".restart",".data")
    fstat_file_name = restart_file_name.replace(".restart",".fstat")

    # open the in- and output files
    data_file = open(data_file_name, 'w')
    fstat_file = open(fstat_file_name, 'w')
    restart_file = open(restart_file_name)
    print("Reading %s" % restart_file_name)
    restart_file = restart_file.readlines()

    #look for number of particles and interactions, get line with first interaction
    i = 0
    while not restart_file[i].startswith('Walls '): i += 1
    num_walls = int(restart_file[i].split()[1])
    index_walls = i + 1
    index_walls_end = index_walls + num_walls

    wall0 = restart_file[index_walls].split()
    torque_index = wall0.index('torque')
    torque_x = 0;
    torque_y = 0;
    torque_z = 0;
    # Read all walls
    for i in range(index_walls, index_walls_end):
        wall = restart_file[i].split()
        wall_id = int(wall[2])
        if wall_id<=wall_id_end and wall_id>=wall_id_begin:
            torque_x += float(wall[torque_index+1])
            torque_y += float(wall[torque_index+2])
            torque_z += float(wall[torque_index+3])
    print("Torque %r %r %r" % (torque_x,torque_y,torque_z))

if __name__ == '__main__':
    main()
