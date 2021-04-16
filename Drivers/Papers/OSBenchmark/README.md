# Comparing open-source DEM frameworks for simulations of common bulk processes

This directory contains the codes used to run the simulations specified in [1] in MercuryDPM.

## Background

The article is based on simulating three different systems:
 - Particles flowing out of a silo
 - A large particle impacting on a bed of small particles
 - A rotating drum filled with a mixture of small and large particles

These three geometries are implemented in MercuryDPM in the following codes:
 - ``OSSilo.cpp``
 - ``OSPenetration.cpp``
 - ``OSDrum.cpp``

The source codes for these three cases in this folder. All three simulations use a shared base class, defined in `MercuryOS.h`, in which the material properties are defined. In addition, there are two self tests that simulate particle pair collisions: ``OSHertzSelfTest.cpp`` and ``OSMindlinSelfTest.cpp``. When you run `make fulltest`, the selftest will be executed and their results will be checked; this will help to detect whether any changes to future versions of MercuryDPM affect the contact model used here.

## Getting started

To run the simulations for this paper, you need to install the developer's version of MercuryDPM (the 'Trunk'). For the data shown in the article, the revision 5520 was used. Here are the installation instructions: https://www.mercurydpm.org/downloads/developers-version-trunk.

For benchmarking, you need to activate OpenMP parallelisation of the code. Therefore, when installling Mercury, you need to run cmake with the option `-DMercury_USE_OpenMP=ON`.

To determine the ideal number of threads to use, run the script ``ompTest.sh``. This will simulate 100 time steps of the silo geometry using different number of threads and return the cpu-time used to run the simulation. Pick the number of threads with the lowest cpu-time and set the variable opt in ``run.sh`` to that number (e.g. ``opt=-omp 12``).

Now execute the script ``run.sh``. This will run all nine simulations used in the paper in sequence.

## References

[1] *Add article reference here*

## Author & Version

Thomas Weinhart, *14 April, 2021*

## License

This project is licensed under the BSD-3 license.