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

# ------------Stereolithography (stl) format output--------------------------


import numpy as np



def make_sphere(pos, R, edges_phi, edges_theta):
    from stl import mesh
    # This function creates an stl model of a sphere
    data = np.zeros((edges_theta - 1) * edges_phi * 2, dtype=mesh.Mesh.dtype)  # Create the data structure for triangles

    x = pos[0]
    y = pos[1]
    z = pos[2]

    k = 0

    for j in range(edges_phi):
        for i in range(edges_theta):
            # v1...v4 - vertices of triangles
            # Triangle 1: v1, v2, v3
            # Triangle 2: v3, v2, v4

            v1 = np.array([x + R * np.sin(np.pi * (i / edges_theta)) * np.cos(2 * np.pi * (j / edges_phi)),
                           y + R * np.sin(np.pi * (i / edges_theta)) * np.sin(2 * np.pi * (j / edges_phi)),
                           z + R * np.cos(np.pi * (i / edges_theta))])
            v2 = np.array([x + R * np.sin(np.pi * ((i + 1) / edges_theta)) * np.cos(2 * np.pi * (j / edges_phi)),
                           y + R * np.sin(np.pi * ((i + 1) / edges_theta)) * np.sin(2 * np.pi * (j / edges_phi)),
                           z + R * np.cos(np.pi * ((i + 1) / edges_theta))])
            v3 = np.array([x + R * np.sin(np.pi * (i / edges_theta)) * np.cos(2 * np.pi * ((j + 1) / edges_phi)),
                           y + R * np.sin(np.pi * (i / edges_theta)) * np.sin(2 * np.pi * ((j + 1) / edges_phi)),
                           z + R * np.cos(np.pi * (i / edges_theta))])
            v4 = np.array([x + R * np.sin(np.pi * ((i + 1) / edges_theta)) * np.cos(2 * np.pi * ((j + 1) / edges_phi)),
                           y + R * np.sin(np.pi * ((i + 1) / edges_theta)) * np.sin(2 * np.pi * ((j + 1) / edges_phi)),
                           z + R * np.cos(np.pi * ((i + 1) / edges_theta))])

            if (not (np.allclose(v1, v3))):
                data['vectors'][k] = np.array([v1, v2, v3])
                k += 1

            if (not (np.allclose(v2, v4))):
                data['vectors'][k] = np.array([v3, v2, v4])
                k += 1

    return data

def save_pebbles_stl(filename, pebbles, N_phi, N_theta):
    from stl import mesh
    data = make_sphere(pebbles[0][:3], pebbles[0][3], N_phi, N_theta)
    for pebble in range(1, len(pebbles)):
        data = np.hstack((data, make_sphere(pebbles[pebble][:3], pebbles[pebble][3], N_phi, N_theta)))
    your_mesh = mesh.Mesh(data, remove_empty_areas=True)
    your_mesh.save(filename)
    print("stl model saved: " + filename)
    return

def save_voxel_grid_stl(filename, vox, bbox, span):
    # A function that saves voxels as an stl
    #   Numbering of vertices v1-v8:
    #        z
    #       |
    #      1----4
    #      |\   |\
    #     |  2----3
    #     5----8  __ y
    #     \    \
    #      6----7
    #      \
    #       x
    N = vox.shape[0]
    from stl import mesh
    data = np.zeros(12 * N ** 3, dtype=mesh.Mesh.dtype)  # Create the data structure for triangles
    l = span / N
    p = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                vox_center = np.array([bbox[0] + l * (1 + 2 * i), bbox[2] + l * (1 + 2 * j), bbox[4] + l * (1 + 2 * k)])
                v1 = vox_center + l * np.array([-1, -1, 1])
                v2 = vox_center + l * np.array([1, -1, 1])
                v3 = vox_center + l * np.array([1, 1, 1])
                v4 = vox_center + l * np.array([-1, 1, 1])
                v5 = vox_center + l * np.array([-1, -1, -1])
                v6 = vox_center + l * np.array([1, -1, -1])
                v7 = vox_center + l * np.array([1, 1, -1])
                v8 = vox_center + l * np.array([-1, 1, -1])

                # Create cube sides
                if (i == 0 and vox[i, j, k]) or (i > 0 and (vox[i, j, k] and not vox[i - 1, j, k])):
                    data['vectors'][p] = np.array([v1, v4, v5])
                    data['vectors'][p + 1] = np.array([v5, v4, v8])
                    p += 2

                if (j == 0 and vox[i, j, k]) or (j > 0 and (vox[i, j, k] and not vox[i, j - 1, k])):
                    data['vectors'][p] = np.array([v6, v2, v1])
                    data['vectors'][p + 1] = np.array([v6, v1, v5])
                    p += 2

                if (k == 0 and vox[i, j, k]) or (k > 0 and (vox[i, j, k] and not vox[i, j, k - 1])):
                    data['vectors'][p] = np.array([v7, v8, v5])
                    data['vectors'][p + 1] = np.array([v7, v5, v6])
                    p += 2

                if (i == N - 1 and vox[i, j, k]) or (i < N - 1 and (vox[i, j, k] and not vox[i + 1, j, k])):
                    data['vectors'][p] = np.array([v3, v6, v2])
                    data['vectors'][p + 1] = np.array([v3, v7, v6])
                    p += 2

                if (j == N - 1 and vox[i, j, k]) or (j < N - 1 and (vox[i, j, k] and not vox[i, j + 1, k])):
                    data['vectors'][p] = np.array([v7, v3, v4])
                    data['vectors'][p + 1] = np.array([v7, v4, v8])
                    p += 2

                if (k == N - 1 and vox[i, j, k]) or (k < N - 1 and (vox[i, j, k] and not vox[i, j, k + 1])):
                    data['vectors'][p] = np.array([v4, v1, v2])
                    data['vectors'][p + 1] = np.array([v4, v2, v3])
                    p += 2

    mesh = mesh.Mesh(data, remove_empty_areas=True)
    mesh.save(filename)
    return

def save_stl_snap(input_mesh, pos, v1, v2, v3, filename):
    print("pos", pos)
    import copy
    ms = copy.deepcopy(input_mesh)
    # Rotate mesh to pd
    e1 = np.array([1,0,0])
    e2 = np.array([0,1,0])
    e3 = np.array([0,0,1])
    Q = np.array([
        [e1@v1, e2@v1, e3@v1],
        [e1@v2, e2@v2, e3@v2],
        [e1@v3, e2@v3, e3@v3]])

    for j in range(len(ms.normals)):
        ms.v0[j]  = Q @ ms.v0[j]
        ms.v1[j]  = Q @ ms.v1[j]
        ms.v2[j]  = Q @ ms.v2[j]
        ms.normals[j]  = Q @ ms.normals[j]

    # Shift Mesh's center of gravity to a specified location
    for j in range(len(ms.normals)):
        ms.v0[j] += pos
        ms.v1[j] += pos
        ms.v2[j] += pos

    ms.save(filename)
    return

def save_stl_sequence(OPT, DATA):
    # Load clump sequence
    filename = OPT['clumpSeqDir'] + 'clump_seq.txt'
    textfile = open(filename, "r")
    content_list = textfile.readlines()

    sequence = np.zeros([len(content_list), 13])
    for i in range(len(content_list)):
        content_list[i] = content_list[i].replace("\n", "")
        word_list = content_list[i].split(' ')
        for j in range (13):
            sequence[i][j] = float(word_list[j])
    DATA['clumpSequence'] = sequence

    # Save stl sequence
    for i in range(len(sequence)):
        pos = sequence[i][1:4]
        v1 = sequence[i][4:7]
        v2 = sequence[i][7:10]
        v3 = sequence[i][10:13]
        print("pos:", pos, "v1:", v1, "v2:", v2, "v3:", v3)
        save_stl_snap(DATA['stlMesh'], pos, v1, v2, v3, './blender/stl_seq/' + 's_' + str(i) + '.stl')


    return OPT, DATA
