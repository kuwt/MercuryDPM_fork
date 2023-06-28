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

# ------------Load data functions----------------------------------------------



def load_pebbles(OPT,DATA):
    import numpy as np
    # This function loads pebbles data in the CLUMP library (csv) format
    filename = OPT['clumpInputDir'] + OPT['clumpFileName']
    textfile = open(filename, "r")
    content_list = textfile.readlines()
    pebbles = np.zeros([len(content_list), 4])
    for i in range(len(content_list)):
        content_list[i] = content_list[i].replace("\n", "")
        word_list = content_list[i].split(',')
        for j in range (4):
            pebbles[i][j] = float(word_list[j])
    DATA['pebbles'] = pebbles
    if (OPT['verbose']): print(str(len(pebbles)) + " pebbles loaded \n")

    return OPT, DATA

def load_mesh(OPT, DATA):
    from stl import mesh
    filename = OPT['stlInputDir'] + OPT['stlFileName']
    stl_mesh = mesh.Mesh.from_file(filename)
    DATA['stlMesh'] = stl_mesh
    return OPT, DATA
