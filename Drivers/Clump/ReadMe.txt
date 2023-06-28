Here you will find the set of examples that demonstrate use cases for multiparticle system of MercuryDPM

Each driver file can be run in two modes

1) "make DriverFile",	"./DriverFile" - manual way of running the examples
2) "./AutoDriverFile" - this executable will erase previous output, re-compile driver file, use data2pvd to create complete paraview output, including energy plots (this way hides all the console output).

Every example is placed in its subdirectory.

Multiparticles (clumps) are generated using MClump tool, located in the MercuryDPM source folder at ./tools/MClump - pls see the readme.txt in that folder.

Please also see a detailed description of the implementation in ./tools/MClump/Doc - the paper 
on the implementation of rigid clumps in MercuryDPM

	
