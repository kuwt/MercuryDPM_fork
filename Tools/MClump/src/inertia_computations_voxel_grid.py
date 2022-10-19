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

# ------------Inertial properties computation based on voxel grid -----------------

import numpy as np
from numba import jit
from src.inertia_computations_pebbles import compute_principal_directions
from src.inertia_computations_pebbles import rotate_to_pd_pebbles_toi

@jit(nopython=True)
def bounding_box(pebbles):
    # returns CUBIC bounding box in the shape [x_min, x_max, y_min, y_max, z_min, z_max]
    # that encloses the clump. This bounding box is used for voxelization. Such an approach
    # is OK if the shapes are not too oblique.

    spans = 0.5 * np.array([np.max(pebbles[:, 0] + pebbles[:, 3]) - np.min(pebbles[:, 0] - pebbles[:, 3]),
                            np.max(pebbles[:, 1] + pebbles[:, 3]) - np.min(pebbles[:, 1] - pebbles[:, 3]),
                            np.max(pebbles[:, 2] + pebbles[:, 3]) - np.min(pebbles[:, 2] - pebbles[:, 3])])
    box_center = 0.5 * np.array([np.max(pebbles[:, 0] + pebbles[:, 3]) + np.min(pebbles[:, 0] - pebbles[:, 3]),
                                 np.max(pebbles[:, 1] + pebbles[:, 3]) + np.min(pebbles[:, 1] - pebbles[:, 3]),
                                 np.max(pebbles[:, 2] + pebbles[:, 3]) + np.min(pebbles[:, 2] - pebbles[:, 3])])
    max_span = np.max(spans)

    return np.array([box_center[0] - max_span,
                     box_center[0] + max_span,
                     box_center[1] - max_span,
                     box_center[1] + max_span,
                     box_center[2] - max_span,
                     box_center[2] + max_span]), max_span

@jit(nopython=True)
def voxel_grid_from_pebbles(pebbles, N):
    # this function returns decomposition of cubic domain with void and material split in binary voxels

    vox = np.zeros(N ** 3).reshape(N, N, N)
    bbox, span = bounding_box(pebbles)

    for i in range(N):
        for j in range(N):
            for k in range(N):
                for pebble in range(len(pebbles)):
                    pos_v = np.array([bbox[0] + 2. * span * (1. / (2. * N) + i / N),
                                      bbox[2] + 2. * span * (1. / (2. * N) + j / N),
                                      bbox[4] + 2. * span * (1. / (2. * N) + k / N)])
                    pos_c = pebbles[pebble, :3]
                    if (np.linalg.norm(pos_c - pos_v) < pebbles[pebble, 3]):
                        vox[i, j, k] = 1.
    return bbox, span, vox

@jit(nopython=True)
def compute_vol_com_voxel_grid(bbox, span, vox):
    # This function defines the absolute position of center of gravity of all voxels
    N = vox.shape[0]
    vol = np.sum(vox) * (2. * span / N)**3

    com = np.array([0.,0.,0.])
    for i in range(N):
        for j in range(N):
            for k in range(N):
                com += vox[i,j,k] * np.array([ bbox[0] + 2.*span *(1./ (2.*N) + i / N),
                                               bbox[2] + 2.*span *(1./ (2.*N) + j / N),
                                               bbox[4] + 2.*span *(1./ (2.*N) + k / N)])
    return vol, com / np.sum(vox)

@jit(nopython=True)
def shift_to_voxel_com_bbox_pebbles(com, pebbles, bbox):
    # Returns the coordinates of pebbles shifted such that the center of mass (COM) is at zero.
    for j in range(len(pebbles)):
        pebbles[j][:3] -= com
    for i in range(3):
        bbox[2 * i] -= com[i]
        bbox[2 * i + 1] -= com[i]

    return pebbles, bbox

@jit(nopython=True)
def clump_toi_voxel_grid(density, com, vox, bbox, span):
    TOI = np.zeros(9).reshape(3,3)
    N = vox.shape[0]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if (vox[i,j,k]):
                    r = np.array([ bbox[0] + 2.*span *(1./ (2.*N) + i / N),
                               bbox[2] + 2.*span *(1./ (2.*N) + j / N),
                               bbox[4] + 2.*span *(1./ (2.*N) + k / N)]) - com
                    X = r[0]
                    Y = r[1]
                    Z = r[2]
                    TOI += np.array([[Y**2 + Z**2, -X*Y, -X*Z],
                                    [-X*Y, X**2 + Z**2, -Y*Z],
                                    [-X*Z, -Y*Z, X**2 + Y**2]])
    return density * (2*span/N)**3 * TOI

def compute_inertia_from_voxel_grid(OPT, DATA):
    # take array of pebbles
    pebbles = DATA['pebbles']
    density = DATA['density']

    # Compute voxel grid
    N = OPT['voxNum']
    bbox, span, vox = voxel_grid_from_pebbles(pebbles, N)

    # Compute mass, and center of mass
    vol, com = compute_vol_com_voxel_grid(bbox, span, vox)

    # Shift pebbles and span to make com an origin
    pebbles, bbox = shift_to_voxel_com_bbox_pebbles(com, pebbles, bbox)

    # Compute tensor of inertia
    toi = clump_toi_voxel_grid(density, com, vox, bbox, span)
    if OPT['verbose']: print("Tensor of inertia of voxels: ", toi)

    # Compute principal directions
    v1, v2, v3 = compute_principal_directions(toi)
    if OPT['verbose']: print("Principal directions (voxel grid): ", v1, v2, v3)

    # Rotate to principal directions
    if OPT['rotateToPD']:
        pebbles, toi, v1, v2, v3 = rotate_to_pd_pebbles_toi(pebbles, toi, v1, v2, v3)
        if OPT['verbose']: print("Rotated principal directions: ", v1, v2, v3)
        if OPT['verbose']:
            print("Rotated toi: ")
            print(toi)

    # Save all the data
    DATA['pebbles'] = pebbles
    DATA['mass'] = vol * DATA['density']
    DATA['pd'] = v1, v2, v3
    DATA['toi'] = toi
    DATA['voxelGrid'] = vox, bbox, span


    return OPT, DATA
