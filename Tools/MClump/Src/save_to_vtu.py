# Copyright (c) 2013-2022, The MercuryDPM Developers Team. All rights reserved.
# For the list of developers, see <http://www.MercuryDPM.org/Team>.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name MercuryDPM nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE

# ------------Unstructured grid (vtu) format output--------------------------


import numpy as np

def save_pebbles_to_vtu(filename, positions, velocities, rad):
    # This function saves pebble information in XML/vtu (unstructured grid) format
    num_part = np.shape(positions)[0]
    num_timesteps = np.shape(positions)[1]

    for timestep in range(num_timesteps):
        pos_string = ""
        vel_string = ""
        rad_string = ""
        for num in range(num_part):
            pos_string = pos_string + "%.4f " % positions[num, timestep, 0]
            pos_string = pos_string + "%.4f " % positions[num, timestep, 1]
            pos_string = pos_string + "%.4f " % positions[num, timestep, 2]
            vel_string = vel_string + "%.4f " % velocities[num, timestep, 0]
            vel_string = vel_string + "%.4f " % velocities[num, timestep, 1]
            vel_string = vel_string + "%.4f " % velocities[num, timestep, 2]
            rad_string = rad_string + "%.4f " % rad[num, timestep]

        pos_string = pos_string + "\n"
        vel_string = vel_string + "\n"
        rad_string = rad_string + "\n"

        textfile = open(filename + '_' + str(timestep) + '.vtu', "w")
        textfile.write('<?xml version="1.0"?>\n')
        textfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        textfile.write(' <UnstructuredGrid>\n')
        textfile.write('  <Piece NumberOfPoints="' + str(num_part) + '" NumberOfCells="0">\n')
        textfile.write('   <Cells>\n')
        textfile.write('    <DataArray type="Int32" name="connectivity" format="ascii">\n')
        textfile.write('       0\n')
        textfile.write('    </DataArray>\n')
        textfile.write('    <DataArray type="Float32" name="offset" format="ascii">\n')
        textfile.write('       0\n')
        textfile.write('    </DataArray>\n')
        textfile.write('    <DataArray type="UInt8" name="types" format="ascii">\n')
        textfile.write('       1\n')
        textfile.write('    </DataArray>\n')
        textfile.write('   </Cells>\n')
        textfile.write('   <Points>\n')
        textfile.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        textfile.write(pos_string)
        textfile.write('</DataArray>\n')
        textfile.write('    </Points>\n')
        textfile.write('    <PointData>\n')
        textfile.write('<DataArray type="Float32" Name="Position" NumberOfComponents="3" format="ascii">\n')
        textfile.write(pos_string)
        textfile.write('</DataArray>\n')
        textfile.write('<DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n')
        textfile.write(vel_string)
        textfile.write('</DataArray>\n')
        textfile.write('<DataArray type="Float32" Name="Radius" NumberOfComponents="1" format="ascii">\n')
        textfile.write(rad_string)
        textfile.write('</DataArray>\n')
        textfile.write('</PointData>\n')
        textfile.write('<CellData/>\n')
        textfile.write('</Piece>\n')
        textfile.write('</UnstructuredGrid>\n')
        textfile.write('</VTKFile>\n')
        textfile.close()

    return




def save_pd_to_vtu(filename, v1, v2, v3):
    # This function saves clump principal directions in XML/vtu (unstructured grid) format
    num_part = 1
    num_timesteps = len(v1)

    for timestep in range(num_timesteps):
        pos_string = ""
        v1_string = ""
        v2_string = ""
        v3_string = ""

        rad_string = ""
        pos_string = pos_string + "%.4f " % 0
        pos_string = pos_string + "%.4f " % 0
        pos_string = pos_string + "%.4f " % 0
        v1_string = v1_string + "%.4f " % v1[timestep, 0]
        v1_string = v1_string + "%.4f " % v1[timestep, 1]
        v1_string = v1_string + "%.4f " % v1[timestep, 2]
        v2_string = v2_string + "%.4f " % v2[timestep, 0]
        v2_string = v2_string + "%.4f " % v2[timestep, 1]
        v2_string = v2_string + "%.4f " % v2[timestep, 2]
        v3_string = v3_string + "%.4f " % v3[timestep, 0]
        v3_string = v3_string + "%.4f " % v3[timestep, 1]
        v3_string = v3_string + "%.4f " % v3[timestep, 2]

        rad_string = rad_string + "%.4f " % 1.0

        pos_string = pos_string + "\n"
        v1_string = v1_string + "\n"
        v2_string = v2_string + "\n"
        v3_string = v3_string + "\n"

        rad_string = rad_string + "\n"

        textfile = open(filename + '_pd_' + str(timestep) + '.vtu', "w")
        textfile.write('<?xml version="1.0"?>\n')
        textfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        textfile.write(' <UnstructuredGrid>\n')
        textfile.write('  <Piece NumberOfPoints="' + str(num_part) + '" NumberOfCells="0">\n')
        textfile.write('   <Cells>\n')
        textfile.write('    <DataArray type="Int32" name="connectivity" format="ascii">\n')
        textfile.write('       0\n')
        textfile.write('    </DataArray>\n')
        textfile.write('    <DataArray type="Float32" name="offset" format="ascii">\n')
        textfile.write('       0\n')
        textfile.write('    </DataArray>\n')
        textfile.write('    <DataArray type="UInt8" name="types" format="ascii">\n')
        textfile.write('       1\n')
        textfile.write('    </DataArray>\n')
        textfile.write('   </Cells>\n')
        textfile.write('   <Points>\n')
        textfile.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        textfile.write(pos_string)
        textfile.write('</DataArray>\n')
        textfile.write('    </Points>\n')
        textfile.write('    <PointData>\n')
        textfile.write('<DataArray type="Float32" Name="Position" NumberOfComponents="3" format="ascii">\n')
        textfile.write(pos_string)
        textfile.write('</DataArray>\n')
        textfile.write('<DataArray type="Float32" Name="v1" NumberOfComponents="3" format="ascii">\n')
        textfile.write(v1_string)
        textfile.write('</DataArray>\n')
        textfile.write('<DataArray type="Float32" Name="v2" NumberOfComponents="3" format="ascii">\n')
        textfile.write(v2_string)
        textfile.write('</DataArray>\n')
        textfile.write('<DataArray type="Float32" Name="v3" NumberOfComponents="3" format="ascii">\n')
        textfile.write(v3_string)
        textfile.write('</DataArray>\n')
        textfile.write('<DataArray type="Float32" Name="Radius" NumberOfComponents="1" format="ascii">\n')
        textfile.write(rad_string)
        textfile.write('</DataArray>\n')
        textfile.write('</PointData>\n')
        textfile.write('<CellData/>\n')
        textfile.write('</Piece>\n')
        textfile.write('</UnstructuredGrid>\n')
        textfile.write('</VTKFile>\n')
        textfile.close()

    return