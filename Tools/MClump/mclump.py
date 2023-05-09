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

#-------------------------------------------------------------------------------
#    MClump----MercuryDPM clump generation tool - main file
#-------------------------------------------------------------------------------


import sys

sys.path.append('./src/')
sys.path.append('./')

try:
    __import__('pip')
except ImportError:
    print("Package manager pip is not installed - impossible to load missing packages")

NumbaInstalled = True
def import_or_install_modules():
    # This function automatically installs required packages if they are missing
    global NumbaInstalled
    try:
        __import__('stl')
    except ImportError:
        print("Package numpy-stl not installed, installing...")
        pip.main(['install', 'numpy-stl'])
    try:
        __import__('numba')
    except ImportError:
        print("Numba not detected (manual installation required), proceed without it")
        NumbaInstalled = False
    return

import_or_install_modules()


from src.output_data import output_clump_data
from src.baseline import baseline
from src.color_lib import colorClass
from src.input_data import load_pebbles
from src.input_data import load_mesh
from src.inertia_computations_pebbles import compute_inertia_from_pebbles
from src.inertia_computations_mesh import compute_inertia_from_mesh
from src.inertia_computations_mixed import compute_inertia_mixed
from src.mesh_to_clump_algorithms import compute_clump_from_mesh

def main():
    # Load baseline parameters - options (OPT) and model data (DATA)
    OPT, DATA = baseline()

    # Set up terminal output
    clr = colorClass()
    clr.set_use_colors(OPT['useColors'])
    clr.def_colors()

    # Hello tag
    print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")
    print(clr.END + clr.BOLD + clr.VIOLET + "     MercuryDPM clump generation tool      ")
    print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")

    # Modification of baseline by the command line arguments (if any)
    if (len(sys.argv) == 3):  # Full format: run -m 2
        if (sys.argv[1] == "-m"):
            OPT['mode'] = int(sys.argv[2])

    # Disable numba usage if enabled, but not installed
    OPT['useNumba'] = OPT['useNumba'] and NumbaInstalled

    if OPT['useNumba']:
        if (OPT['verbose']): print(clr.END + clr.BOLD + clr.GREEN + "Just in time compilation is enabled")
    else:
        if (OPT['verbose']): print(clr.END + clr.BOLD + clr.RED + "Just in time compilation is disabled, voxel grid computations may be slow")



    # Main program
    if OPT['mode']==1:   # 1 - start with the list of pebbles, compute inertia by summation over pebbles
        out = "Mode 1"
        if (OPT['verbose']): out +=": clump input from the list of pebbles, inertial properties via summation over pebbles"
        print(clr.BOLD + clr.GREEN + out)
        # load the input data
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Loading pebble configuration..." + clr.BLUE)
        OPT, DATA = load_pebbles(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

        # Compute mass, center of mass, tensor of inertia, principal directions,
        # shift to center of mass and rotate to principal directions.
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Computing inertial properties..." + clr.BLUE)
        OPT,DATA = compute_inertia_from_pebbles(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

        # compute inertial properties
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Output clump data..." + clr.BLUE)
        OPT, DATA = output_clump_data(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")


    if OPT['mode']==2:   # 2 - start with the list of pebbles, compute inertia by voxelization
        out = "Mode 2"
        if (OPT['verbose']): out += ": clump input from the list of pebbles, inertial properties are via voxel grid"
        print(clr.BOLD + clr.GREEN + out)
        # load the input data
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Loading pebble configuration..." + clr.BLUE)
        OPT, DATA = load_pebbles(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

        # Compute mass, center of mass, tensor of inertia, principal directions,
        # shift to center of mass and rotate to principal directions.

        if OPT['useNumba']:
            from src.inertia_computations_voxel_grid import compute_inertia_from_voxel_grid
            if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Computing inertial properties..." + clr.BLUE)
            OPT, DATA = compute_inertia_from_voxel_grid(OPT, DATA)
            if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

        else:
            from src.inertia_computations_voxel_grid_nonumba import compute_inertia_from_voxel_grid_nonumba
            if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Computing inertial properties..." + clr.BLUE)
            OPT, DATA = compute_inertia_from_voxel_grid_nonumba(OPT, DATA)
            if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")


        # Output clump data
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Output clump data..." + clr.BLUE)
        OPT, DATA = output_clump_data(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

    if OPT['mode']==3:   # 3 - start with the list of pebbles, compute inertia by summation over pebbles
        out = "Mode 3"
        if (OPT['verbose']):
            out += ": clump generation from the stl, inertial properties from stl"
        print(clr.BOLD + clr.GREEN + out)
        if (OPT['verbose']):
            out = clr.BOLD + clr.RED + "WARNING: For better results use external clump generation procedure (e.g CLUMP library)"
        print(clr.BOLD + clr.GREEN + out)

        # load the input data
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Loading stl configuration..." + clr.BLUE)
        OPT, DATA = load_mesh(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

        # Compute mass, center of mass, tensor of inertia, principal directions,
        # shift to center of mass and rotate to principal directions.
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Computing inertial properties..." + clr.BLUE)
        OPT, DATA = compute_inertia_from_mesh(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Computing clump approximation of spheres" + clr.BLUE)
        OPT, DATA = compute_clump_from_mesh(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

        # Output clump data
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Output clump data..." + clr.BLUE)
        OPT, DATA = output_clump_data(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

    if OPT['mode']==4:   # 3 - inertia from stl, external generator of pebbles
        out = "Mode 4"
        if (OPT['verbose']):
            out += ": external clump generation, inertial properties from stl"
        print(clr.BOLD + clr.GREEN + out)

        # load stl mesh
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Loading stl configuration..." + clr.BLUE)
        OPT, DATA = load_mesh(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

        # load pebbles
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Loading pebble configuration..." + clr.BLUE)
        OPT, DATA = load_pebbles(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")

        # Compute mass, center of mass, tensor of inertia, principal directions,
        # shift to center of mass and rotate to principal directions.
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Computing inertial properties..." + clr.BLUE)
        OPT, DATA = compute_inertia_mixed(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")


        # Output clump data
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Output clump data..." + clr.BLUE)
        OPT, DATA = output_clump_data(OPT, DATA)
        if (OPT['verbose']): print(clr.BOLD + clr.YELLOW + "Done")


    print(clr.BOLD + clr.GREEN + "Done!")


if __name__ == "__main__":
    main()

