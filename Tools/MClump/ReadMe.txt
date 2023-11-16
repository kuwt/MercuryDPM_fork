---------------------------------------------
MClump - clump generation tool for MercuryDPM
---------------------------------------------

Usage:

'python mclump.py' - run mclump with default set of parameters
'python mclump.py -m [1..4]' - run mclump with mode specification:

'python mclump.py -m 1' load pebbles, compute inertial properties by summation over pebbles.
'python mclump.py -m 2' load pebbles, compute inertia properties by voxelization
'python mclump.py -m 3' compute inertia properties from stl, load pebbles from file

All the necessary output - mercuryDPM clump properties, vtu unstructured grids, stl surfaces -
will be generated at the specified locations.

Please see ./src/baseline.py for all the tunable parameters and their detailed descriptions

Implementation/requirements:

MClump is written in Python3 and uses numpy/scipy. If numba is installed, some algorithms are
accelerated by just-in-time pre-compilation. MClump uses a third party library numpy-stl - it is
installed automatically by python package manager pip if not installed.

Please see ./doc directory for a detailed up-to-date documentation of MClump functional.




