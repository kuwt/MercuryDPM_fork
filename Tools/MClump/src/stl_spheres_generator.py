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

# This tool generates stl model of multiple spheres, used for test purposes
import numpy as np
def import_or_install_modules():
    # This function automatically installs required packages if they are missing
    try:
        __import__('stl')
    except ImportError:
        print("Package numpy-stl not installed, installing...")
        pip.main(['install', 'numpy-stl'])
    return

import_or_install_modules()



def main():
    from stl import mesh
    from src.save_to_stl import make_sphere
    data = make_sphere([0,0,0], 0.5, 64, 64)
    data = np.hstack((data, make_sphere([-1,0,0], 0.5, 64, 64)))
    data = np.hstack((data, make_sphere([1, 0, 0], 0.5, 64, 64)))

    your_mesh = mesh.Mesh(data, remove_empty_areas=True)
    your_mesh.save('../input/stl/3_spheres_64.stl')
    print("saved")


if __name__ == "__main__":
    main()