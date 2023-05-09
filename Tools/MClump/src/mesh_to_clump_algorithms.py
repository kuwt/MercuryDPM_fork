# Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

# ------------Clump generation algorithms--------------------------------------


import numpy as np
from scipy.spatial import ConvexHull


def bbox_mesh(mesh):
    # Cubic bounding box in the format (x_min, x_max, y_min, y_max, z_min, z_max)
    margin = 0.1
    span_x = 0.5 * (np.max(mesh.v0) - np.min(mesh.v0))
    span_y = 0.5 * (np.max(mesh.v1) - np.min(mesh.v1))
    span_z = 0.5 * (np.max(mesh.v2) - np.min(mesh.v2))

    span = np.max([span_x, span_y, span_z]) + margin

    cent_x = 0.5 * (np.max(mesh.v0) + np.min(mesh.v0))
    cent_y = 0.5 * (np.max(mesh.v1) + np.min(mesh.v1))
    cent_z = 0.5 * (np.max(mesh.v2) + np.min(mesh.v2))

    return np.array([cent_x - span, cent_x + span,
                     cent_y - span, cent_y + span,
                     cent_z - span, cent_z + span]), span

def in_mesh(mesh, point):
    poly = np.vstack((mesh.v0, mesh.v1, mesh.v2))
    hull = ConvexHull(poly)
    new_hull = ConvexHull(np.concatenate((poly, [point])))
    return np.array_equal(new_hull.vertices, hull.vertices)

def regular_placement(mesh, N):
    # Simplest possible algorithm to approximate an stl with regular arrangement of spheres
    bbox, span = bbox_mesh(mesh)
    pebbles = np.zeros([0,4])
    Rad = np.sqrt(3.) * (span / (1. * N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                pos_v = np.array([bbox[0] + 2. * span * (1. / (2. * N) + i / N),
                                  bbox[2] + 2. * span * (1. / (2. * N) + j / N),
                                  bbox[4] + 2. * span * (1. / (2. * N) + k / N)])
                if in_mesh(mesh, pos_v): pebbles = np.vstack((pebbles, np.hstack((pos_v, Rad))))
    return pebbles

def compute_clump_from_mesh(OPT, DATA):
    mesh = DATA['stlMesh']
    N = OPT['regPlacementDefinition']
    if OPT['meshToClumpAlg'] == 1: pebbles = regular_placement(mesh, N)
    DATA['pebbles'] = pebbles
    return OPT, DATA