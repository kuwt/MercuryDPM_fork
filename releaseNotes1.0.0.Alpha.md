# MercuryDPM 1.0.0.Alpha Release notes 	 <img style="float: right;" src="Documentation\Images\mercury-logo-chosen-light.png">

**MercuryDPM 1.0.0.Alpha was released on 1st July 2022.**

**You can find the download instructions here: http://mercurydpm.org -> Downloads -> Release 1.0.0.Alpha**

After 13 years in a self-hosted SVN repository we have moved to git hosted via bitbucket. We only moved the Trunk (now called master) with history kept back to 2018-04-09. This gives you new ways to use and develop the code and we hope this will streamline the process of reintegrating new externally developed features. The old SVN repository will remain read only for a foreseeable future, in case you need an earlier version of MercuryDPM or one of the old releases or development branches. 

The new bitbucket page is found here: https://bitbucket.org/mercurydpm/mercurydpm/ 

---

## REQUIRED DRIVER CHANGES (to update from MercuryDPM 0.11.1)

 - graphviz (dot) is now required to build the documentation 

The following functions changed their name (the line starting with - is the old version, the line starting with + is the new version)

```diff
- setWallsWriteVTK(FileType)
+ setWallsWriteVTK(bool) 
```

```diff
- WallHandler::readTriangleWall(std::string filename, ParticleSpecies* species, Mdouble scaleFactor)
+ WallHandler::readTriangleWall(std::string filename, ParticleSpecies* species, Mdouble scaleFactor, Vec3D shift, Vec3D velocity, Vec3D angularVelocity)
``` 

> BaseParticle is now an abstract class and SphericalParticle is the new class to use in driver codes.
>```diff
>- Baseparticle
>+ SphericalParticle
>```

```diff 
- BaseParticle::getWallInteractionRadius()
+ BaseParticle::getWallInteractionRadius(this)
- BaseParticle::getInteractionRadius()
+ BaseParticle::getInteractionRadius(BaseParticle*)
+ BaseParticle::getMaxInteractionRadius()
```
 
> SuperQuadric class was renamed to SuperQuadricParticle
>```diff
>- SuperQuadric
>+ SuperQuadricParticle
>``` 

```diff
- SuperQuadric::setRadius(Mdouble radius)
+ SuperQuadricParticle::setAxesAndExponents(const Vec3D& axes, const Mdouble& eps1, const Mdouble& eps2)
```

```diff 
- SuperQuadric::computeOverlapAlpha(const LabFixedCoordinates& contactPoint, const LabFixedCoordinates& normal)
+ SuperQuadric::overlapFromContactPoint(const LabFixedCoordinates& contactPoint, const LabFixedCoordinates& normal)
```

```diff 
- DPMBase::broadPhase(BaseParticle* i)
```

```diff
- BaseWall::setVelocityControl(Vec3D forceGoal, Vec3D gainFactor, Vec3D baseVelocity)
+ BaseWall::setForceControl(Vec3D forceGoal, Vec3D gainFactor, Vec3D baseVelocity)
```

```diff 
- setElasticModulusAndRestitutionCoefficient()
+ setEffectiveElasticModulusAndRestitutionCoefficient()
```

```diff
- setShearModulus()
+ setEffectiveShearModulus()
```

```diff
- logger.assert()
+ logger.assert_debug()
``` 
 
```diff 
- BaseClusterInsertionBoundary::setGeometry(Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax)
+ BaseClusterInsertionBoundary::setGeometry(Vec3D posMin, Vec3D posMax)
``` 
 
>Added the possibility for ChuteInsertionBoundary to add a vector of particles, allowing for multiple species in a single InsertionBoundary.
> geometry independent parameters (```radMin```, ```radMax```) moved from ChuteInsertionBoundary to the base class InsertionBoundary.
>```diff 
>- ChuteInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, double radMin, double radMax, double fixedParticleRadius, double inflowVelocity, doubleinflowVelocityVariance)
>+ ChuteInsertionBoundary::set(std::vector<BaseParticle*> particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, double fixedParticleRadius, double inflowVelocity, doubleinflowVelocityVariance)
>+ ChuteInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, double fixedParticleRadius, double inflowVelocity, double inflowVelocityVariance)
>```

>Added the possibility for HopperInsertionBoundary to add a vector of particles, allowing for multiple species in a single InsertionBoundary.
> geometry independent parameters (```radMin```, ```radMax```) moved from HopperInsertionBoundary to the base class InsertionBoundary.
>```diff
>- HopperInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, double yMin, double yMax, double radMin, double radMax, double chuteAngle, double fixedParticleRadius, bool isHopperCentred_, int hopperDim, double hopperAngle, double hopperLength, double hopperExitLength, double hopperHeight, double lift, double fillPercent)
>+ HopperInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, double yMin, double yMax, double chuteAngle, double fixedParticleRadius, bool isHopperCentred_, int hopperDim, double hopperAngle, double hopperLength, double hopperExitLength, double hopperHeight, double lift, double fillPercent)
>+ HopperInsertionBoundary::set(std::vector<BaseParticle*> particleToCopy, unsigned int maxFailed, double yMin, double yMax, double chuteAngle, double fixedParticleRadius, bool isHopperCentred_, int hopperDim, double hopperAngle, double hopperLength, double hopperExitLength, double hopperHeight, double lift, double fillPercent) 
>```
 
> geometry independent parameters (```radMin```, ```radMax```) moved from RandomClusterInsertionBoundary to the base class InsertionBoundary.
>```diff 
>- RandomClusterInsertionBoundary::set(BaseParticle *particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax, Mdouble rMicroParticle)
>+ RandomClusterInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax, Mdouble rMicroParticle)
>- RandomClusterInsertionBoundary::set(BaseParticle &particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax, Mdouble rMicroParticle)
>+ RandomClusterInsertionBoundary::set(BaseParticle& particleToCopy, unsigned int maxFailed, Vec3D posMin,Vec3D posMax, Vec3D velMin, Vec3D velMax, Mdouble rMicroParticle)
>```
>```diff 
>- RandomClusterInsertionBoundary::set(BaseParticle *particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, unsigned int nParticlesPerCluster, Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax)
>+ RandomClusterInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, unsigned int nParticlesPerCluster, Vec3D velMin, Vec3D velMax)
>- RandomClusterInsertionBoundary::set(BaseParticle &particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, unsigned int nParticlesPerCluster, Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax)
>+ RandomClusterInsertionBoundary::set(BaseParticle& particleToCopy, unsigned int maxFailed, Vec3D posMin,Vec3D posMax, unsigned int nParticlesPerCluster, Vec3D velMin, Vec3D velMax)
>```

--- 
 
## KNOWN ISSUES

 - A MPISuperQuadric class which was originally planned to replace MPIParticles when SuperQuadrics are used was added. This class however is not fully functioning at the moment and is replaced by MPIParticle in the whole code for now.
 
 - The Chute and Contraction directory are not up to date and will be resurrected in the future.
 
 - OpenMP version of MercuryCG is not yet implemented.

---

## CHANGES TO DRIVERS

*MercuryMPITests/ForceLawsMPI2Test.cpp*

 - a new test that checks if history parameters like ```slidingSpring``` and ```maxOverlap``` are preserved when a contact jumps from one MPI domain to another

*SelfTests/Boundaries/PSDSelfTest.cpp*

 - a new test for the PSD class
 
*SelfTests/Boundaries/MultiplePSDSelfTest.cpp*

 - a new test for the new insertionBoundary implementation for multiple PSDs in a single InsertionBoundaries
 
*SelfTests/Boundaries/PSDManualInsertionSelfTest.cpp*

 - a new test for the manual insertion routine of the PSD class
 
*SelfTests/Interactions/ConstantRestitutionSelfTest.cpp*

 - a new test for ```BaseSpecies::constantRestitution```
 
*SelfTests/Walls/UnionOfWallsSelfTest.cpp*

 - a new test for BasicUnionOfWalls
 
*UnitTests/ScrewUnitTest.cpp*

 - a new test for the single helix screw wall
 
*MercuryCG/MercuryCG.cpp*

 - added ability to calculate partial fields for mixtures using the command line
 - new option ```-averagebeyonddomain 0``` to ignore particles outside the domain
 
*MercuryCG/CGSelectRegionSelfTest.cpp*

 - a new test for MercuryCG 
 
*SelfTests/Interactions/TorsionFrictionSelfTest.cpp*

 - a new test for FrictionInteraction
 
*SelfTests/Walls/STLRotationSelfTest.cpp*

 - a new test for rotation of STL walls
 
*SelfTests/ComputeVolumeFractionSelfTest.cpp*

 - a new test for Volume fraction computation by CG
 
*MercurySimpleDemos/SelfTestData/HourGlass3DSelfTest.cpp*

 - Corrected SelfTestData
 
*Sinter/SinterBed2SelfTest.cpp*

 - Replaced SelfTestData
 
*Chute/ChuteWithWedge.cpp*

 - added a new driver which considers the flow in a chute past a wedge
 
*SuperQuadricDemos/EllipsoidsBouncingOnWallDemo.cpp*

 - a new driver showing a demo of a non-speherical particle simulation*

*UnitTests/NurbsSurfaceUnitTest.cpp*

*UnitTests/NurbsWallSelfTest.cpp*

*UnitTests/NurbsWallUnitTest.cpp*

 - added new selfTests to test the new NURBS class

*Tutorials/Tutorial10MPI.cpp*

 - added new driver which is a new version of Tutorial10 and shows how to do it in parallel

*Demos/IndustrialMixers*

 - added VerticalMixer, HorizontalMixer and NautaMixer as new demos to show an example of industrial mixing
 
*SimpleDrum/MonodispersedDrum.cpp*

 - added a new directory for drum demos and created a simple demo 

*Template*

 - added a new folder containing a *TemplateDriver.cpp* which can be utilized to write new drivers 

*Papers/MercuryDPM2019/GetDistanceAndNormalForIntersectionOfWalls.cpp*

*Papers/MercuryDPM2019/GetDistanceAndNormalForScrew.cpp*

*Papers/MercuryDPM2019/PlotSilo.cpp*

*Papers/MercuryDPM2019/PlotPolygon.cpp*

*Papers/MercuryDPM2019/PlotNurbs.cpp*

*Papers/MercuryDPM2019/PlotCoil.cpp*

 - added codes to reproduce images for the MercuryDPM2019 paper
 
*SelfTests/Interactions/RollingFrictionSelfTest.cpp*

 - added a new SelfTest for the RollingFriction
 
*SelfTests/Boundaries/CubeDeletionBoundarySelfTest.cpp*

SelfTests/Boundaries/DeletionBoundarySelfTest.cpp*

 - added new SelfTests for the new DeletionBoundary and CubeDeletionBoundary 
 
*Clusters/RandomClusterInsertionBoundaryDemo.cpp*

 - added a demo for a randomized Cluster insertion
 
*UnitTests/FixedClusterInsertionBoundaryUnitTest.cpp*

*UnitTests/RandomClusterInsertionBoundaryUnitTest.cpp*

 - added new UnitTests for the new ClusterInsertionBoundaries 

*SpeedTest/InteractionHandlerSpeedTest.cpp*

*SpeedTest/ParticleHandlerSpeedTest.cpp*

*SpeedTest/SpeedTestWallInteractions.cpp*

 - added new SpeedTests for wall interactions, interaction handling and particle handling
 
*MercurySimpleDemos/FreeCooling2DPeriodicDemo.cpp*

 - added a new demo for Free cooling in a periodic boundary
 
*MercuryMPITests/InsertionBoundaryMPI2Test.cpp*

*MercuryMPITests/LiquidMigrationMPI2Test.cpp*

 - added two new MPI tests for liquid migration and insertion of particles 
 
*SelfTests/Walls/AxisymmetricWallSelfTest.cpp*

 - added a new SelfTest which where an intersection of walls is created in a cylindrical coordinate system
 
*UnitTests/RandomNumbersUnitTest.cpp*

 - random_test was renamed to RandomNumbersUnitTest and moved to the UnitTests directory
 
*UnitTests/LinearViscoelasticUnitTest.cpp*

 - added a new test which tests the LinearViscoelasticNormalSpecies 
 
*MercurySimpleDemos/REV_IsotropicCompression.cpp*

*MercurySimpleDemos/REV_PureShear.cpp*

 - added two simulations of isotropic compression and pure shear in a 3D cuboid box using the StressStrainControlBoundary
 
*SelfTests/Interactions/ChargedBondedInteractionSelfTest.cpp*

 - added a new SelfTest to test the new ```BaseInteraction::isBonded()``` function
 
*Papers/MercuryDPM2019/GetDistanceAndNormalForTriangleWall.cpp*

*Papers/MercuryDPM2019/GetDistanceAndNormalForTriangleWalls.cpp*

 - added two new tests for TriangleWalls, where *GetDistanceAndNormalForTriangleWall.cpp* checks TriangleWalls in a group
 
*TriangleWalls/SingleSphere.cpp*

 - added a new test where a single sphere is rolling over a beam. This was used to test the TriangleWall implementation 
 
*SelfTests/Boundaries/DropletBoundarySelfTest.cpp*

 - added a new SelfTest for the new DropletBoundary
 
*Papers/TunuguntlaBokhoveThornton2014.cpp*

 - added a code for the following paper https://doi.org/10.1017/jfm.2014.223
 
*MercurySimpleDemos/MarbleRun.cpp*

*MercurySimpleDemos/MarbleRunWithSeesaw.cpp*

 - added a basic code to simulate a marble run and a marble run with a seesaw
 
*Papers/OSBenchmark/*

 - added code for the OSBenchmark paper (to be released)
 
*Segregation/SegregationNew.cpp*

*Segregation/Scripts/Analyse.sh*

*Segregation/Matlab/*

 - added a code to simulate segregation on a belt
 
 - added a script to analyse SegregationNew.cpp as well MATLAB files
 
*Segregation/DensitySizeSeg_Belt_WarreJan_25_75.cpp*

 - added belt simulation from a bachelor project  
 
*MercurySimpleDemos/FreeCooling2DinWalls.cpp*

 - added a simple demo where 32^2 particles with the same velocity are placed in a bi-axial box, which makes them collide with the walls and eachother. It is run once with HGrid on and off to test the working (and speedup) of the hgrid*

*SimpleDrum/SquareDrumTest.cpp*

 - added a rolling drum with a squared geometry as demo
 
*SelfTests/Boundaries/PSDManualInsertionSelfTest.cpp*

 - added a new SelfTest for the new manual insertion routine of the PSD class
 
*UnitTests/ArcWallUnitTest.cpp*

 - added a new UnitTest to test the ArcWall feature
 
*Demos/Membrane/MembraneDemo.cpp*

 - added a new demo to demonstrate the membrane feature
 
*SelfTests/Boundaries/NozzleSelfTest.cpp*

 - added a new SelfTest to test the new nozzle feature
 
*Demos/IndustrialMixers/SprayNozzleDemo.cpp*

 - added a new demo to showcase the new liquid drop insertion feature 
 
*SelfTests/Boundaries/DistributionToPSDSelfTest.cpp*

 - added a new SelfTest which creates a discretised normal distribution with the PSD class 
 
*MercurySimpleDemos/LeesEdwardsDemo.cpp*

 - added a demo to showcase the new LeesEdwardsBoundary
 
*MercurySimpleDemos/TimeDependentPeriodicBoundaryTest.cpp*

 - added a new demo to showcase the new TimeDependentPeriodicBoundary feature
 
*MercurySimpleDemos/ShiftingConstantMassFlowMaserBoundarySelfTest.cpp*

 - added a new demo to showcases a shifting ConstantMassFlowMaserBoundary
 
*MercurySimpleDemos/ShiftingMaserBoundarySelfTest.cpp*

 - added a new demo to showcase a shifting SubcriticalMaserBoundary
 
*MercurySimpleDemos/TimeDependentPeriodicBoundary3DSelfTest.cpp*

 - added a new demo to showcase a TimeDependentPeriodicBoundary in 3D
 
--- 
 
## BUG FIXES

#### <u>Major fixes</u>

*Kernel/BaseWall*

 - fixed a bug in ```BaseWall::getInteractionWith``` where particle-wall contact detection was not smooth when transitioning e.g. over a point where two TriangleWall faces touch. Now we keep only the biggest of the overlapping contacts 
 
#### <u>Minor fixes</u>
 
*CMakeModules/MercuryCppFeatureCheck.cmake*

 - changed current check for C++ version from supporting C++11 and C++14 to C++14 only

 - deleted current checks for gcc version 4.9 and 4.8 and replaced it with a check for version 5.0

*CMAKE/FindPython.cmake*

 - fixed Python not working when two versions (Python2 and Python3) where installed. It now prefers Python3 over Python2

 - replaced shebang in .py files to ```usr/bin/env python3```

*Tools/CombineParallelRestartFiles.py*

 - MPI adaption

*Documentation/Pages/CG.dox*

 - fixed typo

*XBalls/xballs.c*

 - removed bug that cause the ```-of``` option to fail

 - fixed that xballs does not exit at keystroke

*Kernel/VTKWriter/SphericalParticleVtkWriter*

 - added the angular velocity to the particle VTU output  

*Kernel/VTKWriter/InteractionVTKWriter*

 - the *InteractionVTKWriter* now correctly specifies the number of components for torque and force

*Kernel/VTKWriter*

 - added a comment line to each *.vtk* file that contains the time value

 - using MPI, *Walls.vtu* files are now only written by processor 0, because all processors create the same file

*Kernel/Walls/TriangleWall*

 - fixed bug in restart information

*Kernel/Walls/Screw*

 - fixed bug in contact detection on screw outer rim (overlap suddenly jumped from 0 to ```thickness_```)

*Kernel/Walls/IntersectionOfWalls*

 - fixed bug in ```convertLimits``` (used in writeVTK)

*Kernel/Walls/AxisymmetricIntersectionOfWalls*

 - fixed bug in ```convertLimits``` (used in writeVTK) 

*Kernel/Walls/BasicIntersectionOfWalls.cc*

 - fixed bug in ```getDistanceAndNormal```

*Kernel/DPMBase *

 - extended the function ```removeOldFiles``` to remove old *.vtu* files in MPI mode as well

 - extended the function ```write``` so it works well in MPI mode as well

 - fixed syntax error when compiling with OMP

 - fixed a bug where VTK output numbering breaks when restarting simulations

*Kernel/Math/Helpers:*

 - fixed bug in ```readOptionalVariable```

*Kernel/Boundaries/StressStrainControlBoundary*

 - fixed bug in ```write``` function

 - fixed bug in ```set``` function, that caused a memory leak

*Kernel/CG/TimeSmoothedCG*

 - fixed bug that caused NaN output

*Kernel/Interactions/FrictionForceInteractions/FrictionInteraction.cc*

 - fixed bug in ```torsionSpring``` implementation: wrong force was computed

*Kernel/Species/Species.h*

 - fixed bug in constructor: BaseSpecies constructor was not called

*Kernel/Math/Helpers.h*

 - fixed bug in ```readOptionalVariable```: the string check didn't work when there was an additional space character before optional variable

*Kernel/Particles/BaseParticle*

 - ```BaseParticle::getInteractionRadius()``` was computed wrong for mixed species. Now there are three functions

   - ```BaseParticle::getMaxInteractionRadius()``` gets the largest possible interaction radius for interactions of this particle with a particle of *any* species

   - ```BaseParticle::getInteractionRadius(BaseParticle*)``` gets the largest possible interaction radius for interactions for a given pair of particles

   - ```BaseParticle::getInteractionRadius(BaseWall*)``` gets the largest possible interaction radius for interactions for a given pair of particle-wall

   - moved computation of the ```interactionRadius``` from the get to the set function. Therefore, two new variables where added:  ```BaseSpecies::interactionRadius``` and ```Species::maxInteractionRadius```

*MATLAB/

 - all scripts are tested and work in *GNU Octave*

*MATLAB/LoadStatistics.m*

 - reworked and moved this file which loads Mercury output files into MATLAB

 - fixed several small bugs

*Kernel/Walls/TriangleWall.cc*

 - fixed bug in ```TriangleWall::getDistanceAndNormal``` which caused the ```edgeBranches``` to be computed wrong

*Kernel/Species/AdhesiveForceSpecies/ReversibleAdhesiveSpecies.cc*

 - fixed a bug in setting of the ```InteractionDistance```. The ```InteractionDistance``` is now only set when there is an ```adhesionStiffness > 0```

*Kernel/Particles/BaseParticle*

 - moved ```ComputeMass()``` here to give each particle type its own ```ComputeMass``` function. This change was necessary to prevent high computational cost from ```dynamic_casts``` of different particle classes (e.g. SuperQuadrics)

*Kernel/Boundaries/PeriodicBoundary*

 - fixed a bug where the HGrid was not updated when particles where shifted due to periodicBoundaries. Now the HGrid is updated every timestep by adding the shift to the ```totalRelativeMaxDisplacement``` of the HGrid

*Tools/MercuryData.h*

 - changed return value of ```MercuryTimeStepIterator<NDIMS>::end()``` in order to fix *data2pvd* tool

*Kernel/Walls/TriangleWall.h*

 - optimized Trianglewall contact detection. In the past there where issues with distinguishing edge and vertex contacts 

*Kernel/Interactions/FrictionForceInteractions/MindlinRollingTorsionInteraction.cc*

 - corrected MindlinRollingTorsionInteraction. ```rollingRelativeVelocity``` is now computed using the radius instead of the diameter

*Drivers\MercuryCG\fstatistics.cpp*

 - *fstatistics* did not work well for the *FiveParticles* example. The problem was that some particles were unduly set as fixed in ```readNextDataFile```. Now, the ```fixed``` property is correctly set, i.e. particles that are ```fixed``` in the restart file (and only those) are ```fixed``` for all time steps

*Kernel/Math/Helpers*

 - ```helpers::readFromCommandLine``` did not read correctly when it's the last argument

*MercuryMPITests/TwoByTwoDomainMPI4Test.cpp*

 - added a new MPI test

*Configuration/Mastermake.cmake*

 - added a maximum test time (timeout) of 80 seconds for SelfTests and UnitTests such that tests will complete even if one hangs

*Kernel/HGrid*

 - fixed a bug in HGrid which was causing a problem on ARM hardware

*Kernel/ParticleHandler*

 - added two missing calls of the function ```DPMBase::handleParticleRemoval``` when particles are removed from the particleHandler

*Documentation/Pages/Installing.dox*

 - fixed a broken link to the documentation's installation instructions

---

## NEW FEATURES

**1. Several old branches, drivers and features were fixed, updated and merged to Trunk.**

**2. Added SelfTestData and InputData to several SelfTests.**

**3. Logger is now well documented and implemented in the whole code to speed up complex outputs.**

 - developers are highly encouraged to use the ```logger()``` function instead of ```std::cout()``` or ```std::cerr()```

 - logger output can now be formatted to a certain width and precision. e.g. ```logger(INFO, "%10.12", normalForce```) will add the normal force with a precision of 10 decimal places and a width of 12.

 - the output of a message with a certain width is now by standard left-aligned. 

 - added a new class object ```doFlush_``` to the logger which prevents it from calling ```std::endl``` when in a logger statement the argument ```Flusher::NOFLUSH``` is called. This is especially useful for developers in write functions where sometimes there are conditional arguments which would else lead to a new invocation of the logger slowing the code.

 - the logger is now responsible for the error handling of the code and has a specific code for certain errors. This is especially useful for coupled simulations where error codes can be passed to other software.



**4. A new PSD class was added to enable the simulation of polydisperse particle systems.**

 - a discretized PSD is now defined by a vector containing radii and probabilities from a separate RadiusAndProbability class. (```std::vector<RadiusAndProbability> particleSizeDistribution_```)

 - an enum class was added to define the TYPE of PSDs. Based on this TYPE, PSDs are converted to the default cumulative number distribution function (CNDF) for the InsertionBoundary class.

 - particles can now be read in from CSV files. Therefore a PSD type object has to be created and the function ```setPSDfromCSV()``` has to be called. ```setPSDfromVector()``` is deprecated and will be removed in future releases.

 - a manual insertion by a fixed volume is now implemented for cases where accuracy for the inserted PSD is of great interest. In this approach particles are inserted class by class from biggest to smallest particles overfilling classes by a probability based approach. This feature is enabled the function ```InsertionBoundary::setManualInsertion()``` to TRUE, default is FALSE.

 - the PSD class can create discretized Normal and Uniform distributions for simulation.

 - the PSD class is now able to compute the first 6 moments of the actual particle distribution.

 - other useful functions are also implemented (```getDx()```, ```ComputeMomenta()```, etc..).

 - additional documentation was added for polydispersity in *Documentation/AdditionalDocumentation*.


**5. The InsertionBoundary was generalized.**

 - particle insertion is solely done by this class while geometry definitions of a boundary are defined in other classes, e.g. CubeInsertionBoundary, HopperInsertionBoundary, etc. Therefore all ```setPSD()```, ```getPSD()```, ```getDistribution()```, ```setDistribution()``` and ```generateParticle()``` functions where moved to the InsertionBoundary class as well as the private members ```psd_```, ```distribution_```, ```radMin_``` and ```radMax_```.

 - for CubeInsertionBoundary the old interface is still active and will generate a uniform distribution on [rMin, rMax] by default. Yet, this will only work for a single species and not a vector of particle species.

 - it is now possible to include several PSDs and Species into a single insertionBoundary.

 - insertionBoundary has a new option to allow insertion without overlap check via a ```CheckParticleForInteraction_``` class member.

 - setters and getters for a PSD and a distribution are now defined. 

 - distribution::Uniform is set as default when no PSD is inserted.

**6. Cleaned and updated documentation for a huge amount of driver codes and kernel features. Moreover, the FAQ for developers has been updated and most of the new features of this list are well documented.**

**7. Walls can now be defined by Non-Uniform Rational B-Splines (NURBS).**

**8. New tutorials are added, old tutorials have been extended and all of them are now well documented.**

**9. A DeletionBoundary and CubeDeletionBoundary is introduced which is able to delete particles being in a specific area at any timestep.**

**10. We added a DropletBoundary which is able to insert liquid droplets. Furthermore, a BoundaryVTKWriter is able to write the droplets into a visualizable .VTK file.**

**11. MercuryDPM can now insert Clusters of particles through the new FixedClusterInsertionBoundary and RandomClusterInsertionBoundary.**

**12. A lot of our C++ code has been reworked and optimized in terms of speed and flexibility.**

 - dynamic_casts where replaced where possible.

 - normal inheritance was preferred over other solutions by adding base classes where possible (this avoids dynamic casting).

 - C++ 14 features where implemented over older C++ standards.

 - output was optimized to flush only when necessary and the whole output is now processed by our own Logger.

 - added an option to run in single precision mode.

 - MPI and OpenMP implementation.

 
**11. Added a MarbleRun simulation to be able to test single particle motion in MercuryDPM.**

**12. A new Membrane class was added to simulate elastic membranes consisting of particles with a certain mass connected by springs.**

 - this class assumes, that the particles are connected in a triangular mesh.

 - a hexagonal unit cell is assumed for a correct computation of the spring constant for the distance of the spring.

 - *Kernel/Walls/MeshTriangle* was added for the membrane class and implements a triangle whose vertex positions are defined by three particles. 


**13. MercuryDPM is now compatible with OpenMP (OMP). An option was added (Mercury_USE_OpenMP, default: OFF) to allow for parallel simulations with OpenMP. This defines the compiler flag MERCURY_USE_OMP, which is used for OMP-specific code.**

 - OMP required the implementation of a reduced algorithms for force computation.

 - ```DPMBase::numberOfOMPthreads_``` stores number threads used when omp is turned on, default: 1.

 - ```DPMBase::setNumberOfOMPthreads``` Warns if OMP is turned of and tests whether OMP is active. If OMP is turned on, values can be set between 1 and ```get_omp_max_threads()```. 


**14 A new Tools directory in the Kernel was created and a csvReader was added to this folder.**

 - the csvReader allows to read in arbitrary data types. For now only reading in the first and the second column for a double type is implemented because it is needed for the PSD class.

 - *.csv* files with a byte-order-mark (BOM) will have the first three bytes of their file skipped to avoid problems (common problem for files created in Windows). All values will be read in as usual.
 
**15. A new contact law for sintering of particles is now part of the Mercury kernel.**
 
**16. Added a new chute class which defines a chute by curvilinear coordinates.**

**17. Added a new TimeDependentPeriodicBoundary class which is able to create a boundary with Lees-Edwards type periodic boundary conditions.**

 - a TimeDependentPeriodicBoundary is like a PeriodicBoundary, but when a particle crosses an edge it is shifted, copied and given a boost.
 
 - See also LeesEdwardsBoundary where particles are only shifted and not copied and given a boost.

 - this type of boundary is useful for studying shear flows.
 
**18. Added a soft stop feature to the Kernel.**

 - after abortion of a code a flag is set which stops the code after the next time step and forces writing of output files.
 
**19. Several smaller features added:**

*Kernel/Walls/Screw*

 - generalised Screw; it can now create both single and double helices
 
*Kernel/Walls/BasicUnionOfWalls*

 - added new Wall type to define a union of two walls. The difference between two walls and a union is that a contact that flips from one wall to the other is preserved for a union, whereas it is considered a new contact for two separate walls

*Kernel/Walls/BaseWall*

 - added a function ```setVelocityControl``` that can be used to apply a certain pressure on a wall; the control is very simple and should be improved

 - added a function ```addRenderedWall``` that allows you to specify a wall to be rendered in place of the actual wall. The renderer applies the position and orientation of the actual wall to the rendered wall, so the rendered wall moves like the original wall

 - added function ```addParticlesAtWall```, which allows you to quickly add particles touching a wall (this was previously done only in Driver codes)

*Kernel/Species/BaseSpecies and Interaction/BaseInteraction*

 - added a property ```constantRestitution``` (FALSE by default). If TRUE, the normal contact law is modified to produce a restitution coefficient and collision time that is independent of mass (by making stiffness and dissipation proportional to mass). This only works right now for the LinearPlastic contact law

*Kernel/Interaction/BaseInteraction*

 - generalised ```BaseInteraction::getEffectiveRadius``` and ```BaseInteraction::getEffectiveMass``` to use invMass and curvature instead of mass and radius; now the computation is the same for any particle or wall type

 - added new function ```BaseInteraction::isBonded()``` which returns TRUE if the interaction is an 'internal' bond. these internal forces are not written to fstat and thus don't contribute to fstatistics

*Kernel/BaseInteractable*

 - added functions ```getInvMass``` and ```getCurvature```, also for walls

*Kernel/Species/LinearViscoelasticNormalSpecies*

*Kernel/Species/LinearPlasticViscoelasticNormalSpecies*

 - added a function ```setRestitutionCoefficient``` that sets the restitution coefficient while keeping stiffness constant

 - added a function ```computeBondNumberMax``` to easily get the max Bond number that this contact law can produce

 - added a new option ```doConstantUnloadingStiffness()``` to behave like Walton-Braun

*Kernel/MercuryTime.h*

 - added functions ```tic``` and ```toc```, which allows you to quickly profile a code's speed

 - added a new function ```getWallTime()``` to get the main simulation time (time of a clock on the wall) and ```getCPUTime()``` to get pure CPU computation time of a simulation

*Kernel/Logger*

 - added colour to logger errors in MPI

*Kernel/Math/Matrix*

 - added new function diag, returning the diagonal elements

*Kernel/Math/Vector*

 - added new functions ```multiplyElementwise``` and ```signedSquare```

*Kernel/Math/ExtendedMath.h*

 - a new constant ```constants::degree``` has been defined as a shortcut for pi/180

*Kernel/Math/RNG*

 - moved seed function out of the loop, so it's only called once

*Kernel/Math/Helpers:*

 - added new function ```readFromCommandLine```, to make it easy to interpret command line arguments

 - added new function ```writeCommandLineToFile```, so you can store command line arguments

 - added new function ```linspace(double, double, int)``` which works the same way as in MATLAB

 - added new function ```removeFromCommandline(int& argc, char* argv[], std::string varName, int nArgs)```. This function  does effectively remove specified commandline arguments from argc and argv

*Kernel/GeneralDefine*

 - redefined ```NUMBER_OF_PROCESSORS``` and ```PROCESSOR_ID``` such that they are also set (as 1 and 0) for serial codes; thus, you can write one code that can be run in serial or parallel, without the use of preprocessor directives

*Kernel/File*

 - added function ```writeFirstAndLastTimeStep``` that sets the saveCount to infinity

*Kernel/DPMBase*

 - added parameter ```readSpeciesFromDataFile```;  if set to TRUE, then the info column in the data file will be interpreted as the species during restarting

 - added function ```fillDomainWithParticles``` created an initially filled geometry (allows a quick check if the geometry is implemented correctly)

 - added function ```setInteractionsWriteVTK(bool)```

 - added function ```get/setReadInteractions()```

 - added function ```DPMBase::getMomentum()``` and ```DPMBase::getAngularMomentum()```

 - added function ```splitDomain(DomainSplit)```, which picks a sensible domain composition for the number of processors (e.g. ```splitDomain(XY)``` splits domain into a 6x5x1 grid for 30 processors)

 - added options to ```readRestartFile()``` which allows to read a restart file without interactions and/or particles

 - Merucry now records the time when the simulation starts

 - added time measurement of simulations. ```clock_.tic()``` starts measuring and ```clock_toc()``` ends measuring time

*Kernel/CG/CGHandler*

 - CGHandler now allows setting an initial file counter (to start reading from e.g. **.data.100*)

*Kernel/CG/BaseCG*

 - added functions ```set[X,Y,Z]``` and ```Grid(min,max,h)```, that lets you quickly define a grid of given mesh size

*Kernel/Boundaries/LeesEdwardsBoundary*

 - added function ```updateBoundaries``` that lets you update the position of a LeesEdwardsBoundary, instead of completely resetting it

*Kernel/Boundaries/InsertionBoundary*

 - added function ```insertParticles```, that lets you insert particles once, without adding the boundary to the handler

 - added functions ```setInitialVolume(Mdouble initialVolume)```, ```getInitialVolume()``` and ```setVolumeFlowRate_(Mdouble volumeFlowRate)```, ```getVolumeFlowRate_()``` that let's you set an initial volume and a constant flow rate for insertion

 - added function ```setVariableVolumeFlowRate```, that lets you define a variable flow rate

 - added ```psd``` and ```distribution``` to ```InsertionBoundary::write()```

*CMakeModules/MercuryXballs.cmake*

 - turning off xballs support for Windows by default 

*Kernel/Walls/TriangleWall.h*

 - The position of a TriangleWall can now be set explicitly

*Kernel/Particles/SuperQuadricParticle*

 - the radius of a SuperQuadric is now used to store the bounding radius. This makes contact detection more efficient

 - setting the ```BoundingRadius``` was reworked to be implmented in a more general way

 - superquadrics without contact were added to the interactionHandler, now superquadrics without contacts are removed from the interactionHandler 

*Kernel/Particles/BaseParticle*

 - added ```BaseParticle::getMomentum()``` and ```BaseParticle::getAngularMomentum()``` to automate some computations

*Kernel/BaseObject*

 - added ```groupID``` which makes it easier to apply an action to a specific group in a handler, visualise a specific group, coarse-grain a specific group, etc.

*Kernel/Particles/LiquidFilmParticle*

 - added ```liquidFilmVolume``` and ```liquidBridgeVolume``` to vtk output

*Kernel/Boundaries/DeletionBoundary*

 - added feature to track the outflowing particles

 - added copy constructor

*Kernel/MpiDataClass:*

 - added new function ```getMPISum()``` which sums the values over all processors

 - MPIParticle can now be changed from SphericalParticle to LiquidFilmParticle with a simple one-line change

*Kernel/InteractionHandler*

 - added a new function ```getLiquidBridgeVolume()```

*Tools/CombineParallelDataFiles.cpp*

 - added to create single data file from all data files created by each process when run in parallel

*Tools/CombineParallelOutFiles.py*

 - added python script to merge .out files

*Tools/CombineParallelOutFiles.py*

 - added python script to merge .out files

*Kernel/Mercury3D*

 - added distance check in ```Mercury3D:hGridFindContactsWithTargetCell()``` which greatly speeds up simulations

*Kernel/GeneralDefine*

 - added ```NaN```, ```inf```, ```intMax``` and ```unsignedMax``` to namespace constants

 - added gas constant ```R``` to namespace constants

*Documentation/AdditionalDocumentation/LoggerImplementationDetail.md*

 - added additional documentation for the logger

*XBalls/xballs.c*

 - updated XBalls

*Kernel/Particles/HeatFluidCoupledParticle*

 - added a new particle type that stores both temperature/heat capacity and liquid content which is adapted for the CFD-DEM studies

*Kernel/Interactions/NormalForceInteractions/SinterLinInteraction.cc*

*Kernel/Interactions/NormalForceInteractions/SinterLinInteraction.h*

*Kernel/Species/NormalForceSpecies/SinterLinNormalSpecies.cc*

*Kernel/Species/NormalForceSpecies/SinterLinNormalSpecies.h*

*Kernel/Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h*

 - added new species and interaction for the new sintering contact model

*Configuration/Docker/dockerfile*

 - added a dockerfile to the configuration folder. This dockerfile is used to get Pipelines on Bitbucket running, but can also be used for other purposes
