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

# ------------Output data functions---------------------------


from src.save_to_vtu import save_pebbles_to_vtu
from src.save_to_vtu import save_pd_to_vtu
from src.save_to_stl import save_voxel_grid_stl
from src.save_to_stl import save_pebbles_stl

import numpy as np
from src.inertia_computations_pebbles import rotate_to_pd_pebbles

def paraview_output_pebbles_animated(OPT, DATA):
    import os
    pebbles = DATA['pebbles']
    num_timesteps = OPT['numTimesteps']

    filename = './clumps/' + OPT['clumpName']
    if not os.path.exists(filename):
        os.mkdir(filename)

    filename += '/' + OPT['paraviewOutputDir']
    if not os.path.exists(filename):
        os.mkdir(filename)

    filename += '/' + OPT['vtuFileName']
    if not os.path.exists(filename):
        os.mkdir(filename)

    positions = np.zeros([len(pebbles), num_timesteps, 3])
    velocities = np.zeros([len(pebbles), num_timesteps, 3])
    v1t = np.zeros([num_timesteps, 3])
    v2t = np.zeros([num_timesteps, 3])
    v3t = np.zeros([num_timesteps, 3])

    radii = np.zeros([len(pebbles), num_timesteps])

    for k in range(num_timesteps):

        v1 = np.array([np.cos(k * (2. * np.pi / (num_timesteps + 0.))),  -np.sin(k * (2. * np.pi / (num_timesteps + 0.))), 0.0])
        v2 = np.array([np.sin(k * (2. * np.pi / (num_timesteps + 0.))), np.cos(k * (2. * np.pi / (num_timesteps + 0.))), 0.0])
        v3 = np.array([0.0, 0.0, 1.0])

        v1c = np.array([np.cos(-k * (2. * np.pi / (num_timesteps + 0.))), -np.sin(-k * (2. * np.pi / (num_timesteps + 0.))), 0.0])
        v2c = np.array([np.sin(-k * (2. * np.pi / (num_timesteps + 0.))), np.cos(-k * (2. * np.pi / (num_timesteps + 0.))), 0.0])
        v3c = np.array([0.0, 0.0, 1.0])

        r_pebbles = rotate_to_pd_pebbles(pebbles, v1, v2, v3)
        v1t[k][:] = v1c
        v2t[k][:] = v2c
        v3t[k][:] = v3c

        for i in range(len(pebbles)):
            positions[i][k][:] = r_pebbles[i][:3]
            radii[i][k] = r_pebbles[i][3]

    save_pebbles_to_vtu(filename, positions, velocities, radii)
    save_pd_to_vtu(filename, v1t, v2t, v3t)

    return

def paraview_output_pebbles_static(OPT, DATA):
    import os
    pebbles = DATA['pebbles']

    filename = './clumps/' + OPT['clumpName']
    if not os.path.exists(filename):
        os.mkdir(filename)

    filename += '/' + OPT['paraviewOutputDir']
    if not os.path.exists(filename):
        os.mkdir(filename)

    filename += '/' + OPT['vtuFileName']
    if not os.path.exists(filename):
        os.mkdir(filename)

    if os.path.exists(OPT['paraviewOutputDir']):  # Remove vtu files from the previous run
        os.chdir(OPT['paraviewOutputDir'])
        os.remove('*.vtu')

    positions = np.zeros([len(pebbles), 1, 3])
    velocities = np.zeros([len(pebbles), 1, 3])
    radii = np.zeros([len(pebbles), 1])

    for i in range(len(pebbles)):
        positions[i][0][:] = pebbles[i][:3]
        radii[i][0] = pebbles[i][3]

    save_pebbles_to_vtu(filename, positions, velocities, radii)
    return


def save_clump_mdpm(OPT, DATA):
    import os

    path = './clumps/' + OPT['clumpName']
    if not os.path.exists(path):
        os.mkdir(path)

    pathI = path + '/' + OPT['inertiaOutputDir'] + '/'
    if not os.path.exists(pathI):
        os.mkdir(pathI)

    pathC = path + '/' + OPT['clumpOutputDir'] + '/'
    if not os.path.exists(pathC):
        os.mkdir(pathC)

    np.savetxt(pathI + OPT['inertiaFileName'], DATA['toi'], delimiter=',')
    np.savetxt(pathI + OPT['pdFileName'], DATA['pd'], delimiter=',')
    np.savetxt(pathI + OPT['massFileName'], [DATA['mass']], delimiter=',')
    np.savetxt(pathC + OPT['outClumpFileName'], DATA['pebbles'], delimiter=',')
    return OPT, DATA


def output_clump_data(OPT, DATA):
    # Output clump data for MDPM + additional output (vtu, stl)
    import os

    # MDPM output
    OPT, DATA = save_clump_mdpm(OPT, DATA)

    # Paraview output
    if OPT['paraviewOutput']:
        if OPT['verbose']: print("Output vtu data")
        if (OPT['animate']): paraview_output_pebbles_animated(OPT, DATA)
        if not (OPT['animate']): paraview_output_pebbles_static(OPT, DATA)

    # Pebble stl output
    if OPT['PebbleStlOutput']:
        if OPT['verbose']: print("Output pebble stl data")
        path = './clumps/' + OPT['clumpName']
        if not os.path.exists(path):
            os.mkdir(path)

        path += '/' + OPT['stlOutputDir'] + '/'
        if not os.path.exists(path):
            os.mkdir(path)

        filename = path + OPT['outStlFileName']
        pebbles = DATA['pebbles']
        N_theta = OPT['stlThetaRes']
        N_phi = OPT['stlPhiRes']
        save_pebbles_stl(filename, pebbles, N_theta, N_phi)

    # Voxel grid stl output
    if OPT['mode'] == 2 and OPT['VoxelStlOutput']:
        if OPT['verbose']: print("Output voxel grid stl data")
        path = './clumps/' + OPT['clumpName']
        if not os.path.exists(path):
            os.mkdir(path)

        path += '/' + OPT['voxelStlOutputDir'] + '/'
        if not os.path.exists(path):
            os.mkdir(path)

        filename = path + OPT['voxStlFileName']
        vox, bbox, span = DATA['voxelGrid']
        save_voxel_grid_stl(filename, vox, bbox, span)

    return OPT, DATA

