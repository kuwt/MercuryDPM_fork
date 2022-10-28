Here you will find the set of examples that demonstrate use cases for multiparticle system of MercuryDPM

Each driver file can be run in two modes

1) "make driver_file",	"./driver file" - manual way of running the examples
2) "./auto_driver_file" - this executable will erase previous output, re-compile driver file, use data2pvd to create complete paraview output, including energy plots (this way hides all the output).

Example 1: Single.cpp - shows energy equipartition in a system of a single nonspherical particle confined in the elastic box.

Example 2: TBar.cpp - shows Dzhanibekov effect on the example of a T-bar rotating around its intermediate axis

Example 3: BulkTs.cpp - a pile of T-bars in a gravity field in presence of dissipation

Example 4: Gomboc.cpp - a rolly - polly out of simply-commected shape of constant density

Multiparticles (clumps) are generated using MClump tool, located in the MercuryDPM source folder at ./tools/MClump - pls see the readme.txt in that folder.

Please also see a detailed description of the implementation in ./tools/MClump/Doc - updated concurrently with the development

	
