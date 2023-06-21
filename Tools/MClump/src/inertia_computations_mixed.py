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
	
# ------------Inertial properties computation based on stl mesh-----------------

import numpy as np
from src.inertia_computations_pebbles import compute_principal_directions


def compute_mass_com_mesh(mesh, density):
    # Computes center of mass of the body bound by the triangulated surface "mesh"
    # Returns mass, coordinates of the center of mass
    O = np.zeros(3)
    com = 0
    vol = 0
    for j in range(len(mesh.normals)):
        p1 = mesh.v0[j]
        p2 = mesh.v1[j]
        p3 = mesh.v2[j]
        vd = np.vstack((np.hstack((O,1)),
                           np.hstack((p1,1)),
                           np.hstack((p2,1)),
                           np.hstack((p3,1))))    
        V_j = -(1./6.) * np.linalg.det(vd)
        c_j = (O + p1 + p2 + p3)/4.
        com += V_j * c_j
        vol += V_j
    com /= vol
    return density*vol, com


def shift_to_center_of_mass_mesh_pebbles(mesh, pebbles, density):
    # Returns the coordinates of mesh and pebbles shifted such that the center of mass is at zero.
    mass, com = compute_mass_com_mesh(mesh, density)
    for j in range(len(mesh.normals)):
        mesh.v0[j] -= com
        mesh.v1[j] -= com
        mesh.v2[j] -= com

    for j in range(len(pebbles)):
        pebbles[j][:3] -= com

    return mass, mesh, pebbles


def clump_toi_from_mesh(mesh, density):
    # This function computes the tensor of inertia of a body bound by a triangulated surface
    O = np.zeros([3])
    a = 0
    b = 0
    c = 0
    a_p = 0
    b_p = 0
    c_p = 0
    for j in range(len(mesh.normals)):
        x_1 = O[0]
        y_1 = O[1]
        z_1 = O[2]
        x_2 = mesh.v0[j][0]
        y_2 = mesh.v0[j][1]
        z_2 = mesh.v0[j][2]
        x_3 = mesh.v1[j][0]
        y_3 = mesh.v1[j][1]
        z_3 = mesh.v1[j][2]
        x_4 = mesh.v2[j][0]
        y_4 = mesh.v2[j][1]
        z_4 = mesh.v2[j][2]
        
        vd = np.array(  [  [x_1, y_1, z_1, 1],
                           [x_2, y_2, z_2, 1],
                           [x_3, y_3, z_3, 1],
                           [x_4, y_4, z_4, 1]   ]) 
        V_j = - (1./6.) * np.linalg.det(vd)
        
        a += V_j * ( y_1**2 + y_1*y_2 + y_2**2 + y_1*y_3 + y_2*y_3 +
                y_3**2 + y_1*y_4 + y_2*y_4 + y_3*y_4 + y_4**2 + 
                z_1**2 + z_1*z_2 + z_2**2 + z_1*z_3 + z_2*z_3 + 
                z_3**2 + z_1*z_4 + z_2*z_4 + z_3*z_4 + z_4**2 ) / 10.
    
        b += V_j * ( x_1**2 + x_1*x_2 + x_2**2 + x_1*x_3 + x_2*x_3 + 
                x_3**2 + x_1*x_4 + x_2*x_4 + x_3*x_4 + x_4**2 + 
                z_1**2 + z_1*z_2 + z_2**2 + z_1*z_3 + z_2*z_3 + 
                z_3**2 + z_1*z_4 + z_2*z_4 + z_3*z_4 + z_4**2 ) / 10.
    
        c += V_j * ( x_1**2 + x_1*x_2 + x_2**2 + x_1*x_3 + x_2*x_3 + 
                x_3**2 + x_1*x_4 + x_2*x_4 + x_3*x_4 +  x_4**2 + 
                y_1**2 + y_1*y_2 + y_2**2 + y_1*y_3 + y_2*y_3 + 
                y_3**2 + y_1*y_4 + y_2*y_4 + y_3*y_4 + y_4**2 ) / 10.
    
        a_p += V_j * ( 2*y_1*z_1 + y_2*z_1 + y_3*z_1 + y_4*z_1 + y_1*z_2 + 
                  2*y_2*z_2 + y_3*z_2 + y_4*z_2 + y_1*z_3 + y_2*z_3 + 
                  2*y_3*z_3 + y_4*z_3 + y_1*z_4 + y_2*z_4 + y_3*z_4 + 
                  2*y_4*z_4) / 20.
    
        b_p += V_j * ( 2*x_1*z_1 + x_2*z_1 + x_3*z_1 + x_4*z_1 + x_1*z_2 + 
                  2*x_2*z_2 + x_3*z_2 + x_4*z_2 + x_1*z_3 + x_2*z_3 + 
                  2*x_3*z_3 + x_4*z_3 + x_1*z_4 + x_2*z_4 + x_3*z_4 + 
                  2*x_4*z_4) / 20.
    
        c_p += V_j * ( 2*x_1*y_1 + x_2*y_1 + x_3*y_1 + x_4*y_1 + x_1*y_2 + 
                  2*x_2*y_2 + x_3*y_2 + x_4*y_2 + x_1*y_3 + x_2*y_3 + 
                  2*x_3*y_3 + x_4*y_3 + x_1*y_4 + x_2*y_4 + x_3*y_4 + 
                  2*x_4*y_4) / 20.
    
    return density * np.array([[a, -b_p, -c_p],[ -b_p, b, -a_p],[ -c_p, -a_p, c]])


def rotate_to_principal_directions_mesh_toi_pebbles(mesh, toi, pebbles, v1, v2, v3):
    # This function rotates the centered mesh to match specified principal directions v1, v2, v3 with Cartesian axes
    e1 = np.array([1,0,0])
    e2 = np.array([0,1,0])
    e3 = np.array([0,0,1])
    Q = np.array([
    [e1@v1, e2@v1, e3@v1],
    [e1@v2, e2@v2, e3@v2],
    [e1@v3, e2@v3, e3@v3]])
    
    for j in range(len(mesh.normals)):
        mesh.v0[j]  = Q @ mesh.v0[j]
        mesh.v1[j]  = Q @ mesh.v1[j]
        mesh.v2[j]  = Q @ mesh.v2[j]
        mesh.normals[j]  = Q @ mesh.normals[j]
    toi = np.transpose(Q) @ toi @ Q

    r_pebbles = np.zeros(np.shape(pebbles))
    for j in range(len(pebbles)):
        r_pebbles[j][:3] = Q @ pebbles[j][:3]
        r_pebbles[j][3] = pebbles[j][3]

    return mesh, toi, r_pebbles, e1, e2, e3


# Extra rotation of the mesh
def rotate_to_pd(mesh, v1, v2, v3):
    # This function rotates the centered mesh to match specified principal directions v1, v2, v3 with Cartesian axes
    e1 = np.array([1,0,0])
    e2 = np.array([0,1,0])
    e3 = np.array([0,0,1])
    Q = np.array([
        [e1@v1, e2@v1, e3@v1],
        [e1@v2, e2@v2, e3@v2],
        [e1@v3, e2@v3, e3@v3]])

    for j in range(len(mesh.normals)):
        mesh.v0[j]  = Q @ mesh.v0[j]
        mesh.v1[j]  = Q @ mesh.v1[j]
        mesh.v2[j]  = Q @ mesh.v2[j]
        mesh.normals[j]  = Q @ mesh.normals[j]

    return mesh


def compute_inertia_mixed(OPT, DATA):
    # take arrays of mesh and pebbles
    mesh = DATA['stlMesh']
    density = DATA['density']
    pebbles = DATA['pebbles']
    # Compute mass, shift to center of mass
    mass, mesh, pebbles = shift_to_center_of_mass_mesh_pebbles(mesh, pebbles, density)
    if OPT['verbose']: print("Total mass of stl: ", mass)

    # Compute tensor of inertia
    toi = clump_toi_from_mesh(mesh, density)
    if OPT['verbose']: print("Tensor of inertia of stl: ", toi)

    # Compute principal directions
    v1, v2, v3 = compute_principal_directions(toi)
    if OPT['verbose']: print("Principal directions: ", v1, v2, v3)

    # Rotate to principal directions
    if OPT['rotateToPD']:
        mesh, toi, pebbles, v1, v2, v3 = rotate_to_principal_directions_mesh_toi_pebbles(mesh, toi, pebbles, v1, v2, v3)
        if OPT['verbose']: print("Rotated principal directions: ", v1, v2, v3)
        if OPT['verbose']: print("Rotated toi: ", toi)

    #mesh = rotate_to_pd(mesh, np.array([1,0,0]), np.array([0,0,-1]), np.array([0,1,0]))


    # Save all the data
    DATA['stlMesh'] = mesh
    DATA['pebbles'] = pebbles
    DATA['mass'] = mass
    DATA['pd'] = v1, v2, v3
    DATA['toi'] = toi

    return OPT, DATA
