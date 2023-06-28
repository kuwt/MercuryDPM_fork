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

# ------------Inertial properties computation based on the list of pebbles-----

import numpy as np


def ComputeMassCOMPebbles(pebbles, density):
    # Computes mass/center of mass (COM) of the assembly of pebbles
    com = np.zeros(3)
    mass = 0
    for j in range(len(pebbles)):
        c_j = pebbles[j][:3]
        r = pebbles[j][3]
        m_j = (4. / 3.) * np.pi * r ** 3 * density
        com += m_j * c_j
        mass += m_j
    com /= mass
    return mass, com


def ShiftToCOMPebbles(pebbles, density):
    # Returns the coordinates of pebbles shifted such that the center of mass (COM) is at zero.
    mass, com = ComputeMassCOMPebbles(pebbles, density)
    for j in range(len(pebbles)):
        pebbles[j][:3] -= com
    return mass, pebbles


def ComputeTOIPebbles(pebbles, density):
    # Tensor of inertia of non-overlapping spherical particles (origin of CS is in center of mass)
    toi = np.zeros(9).reshape(3, 3)
    for j in range(len(pebbles)):
        X = pebbles[j][0]
        Y = pebbles[j][1]
        Z = pebbles[j][2]
        r = pebbles[j][3]
        m = density * (4. / 3.) * np.pi * r ** 3

        toi += m * np.array([[0.4 * r ** 2 + Y ** 2 + Z ** 2, -X * Y, -X * Z],
                             [-X * Y, 0.4 * r ** 2 + X ** 2 + Z ** 2, -Y * Z],
                             [-X * Z, -Y * Z, 0.4 * r ** 2 + X ** 2 + Y ** 2]])
    return toi


def ComputePrincipalDirections(toi):
    # For a given tensor of inertia "toi" returns normalized principal directions
    w, v = np.linalg.eig(toi)
    v1 = v[0] / np.linalg.norm(v[0])
    v2 = v[1] / np.linalg.norm(v[1])
    v3 = v[2] / np.linalg.norm(v[2])

    # Make sure found PDs form the RIGHT basis
    if np.allclose(np.cross(v1, v2), -v3): v3 = -v3
    return v1, v2, v3

def RotateToPDPebbles(pebbles, v1, v2, v3):
    # This function rotates the centered pebbles to the specified principal directions v1, v2, v3.
    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])
    Q = np.array([
        [e1 @ v1, e2 @ v1, e3 @ v1],
        [e1 @ v2, e2 @ v2, e3 @ v2],
        [e1 @ v3, e2 @ v3, e3 @ v3]])

    r_pebbles = np.zeros(np.shape(pebbles))
    for j in range(len(pebbles)):
        r_pebbles[j][:3] = Q @ pebbles[j][:3]
        r_pebbles[j][3] = pebbles[j][3]

    return r_pebbles



def RotateToPDPebblesTOI(pebbles, toi, v1, v2, v3):
    # This function rotates the centered pebbles to the specified principal directions v1, v2, v3.
    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])
    Q = np.array([
        [e1 @ v1, e2 @ v1, e3 @ v1],
        [e1 @ v2, e2 @ v2, e3 @ v2],
        [e1 @ v3, e2 @ v3, e3 @ v3]])

    for j in range(len(pebbles)):
        pebbles[j][:3] = Q @ pebbles[j][:3]

    toi = np.transpose(Q) @ toi @ Q
    return pebbles, toi, e1, e2, e3


def ComputeInertiaFromPebbles(OPT, DATA):
    # take array of pebbles
    pebbles = DATA['pebbles']
    density = DATA['density']
    # Compute mass, shift to center of mass
    mass, pebbles = ShiftToCOMPebbles(pebbles, density)
    if OPT['verbose']: print("Total mass of pebbles: ", mass)

    # Compute tensor of inertia
    toi = ComputeTOIPebbles(pebbles, density)
    if OPT['verbose']: print("Tensor of inertia of pebbles: ", toi)

    # Compute principal directions
    v1, v2, v3 = ComputePrincipalDirections(toi)
    if OPT['verbose']: print("Principal directions: ", v1, v2, v3)

    # Rotate to principal directions
    if OPT['rotateToPD']:
        pebbles, toi, v1, v2, v3 = RotateToPDPebblesTOI(pebbles, toi, v1, v2, v3)
        if OPT['verbose']: print("Rotated principal directions: ", v1, v2, v3)
        if OPT['verbose']: print("Rotated toi: ", toi)

    # Save all the data
    DATA['pebbles'] = pebbles
    DATA['mass'] = mass
    DATA['pd'] = v1, v2, v3
    DATA['toi'] = toi


    return OPT, DATA
