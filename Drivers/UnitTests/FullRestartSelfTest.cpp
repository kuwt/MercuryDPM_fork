//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "DPMBase.h"
#include "Walls/Coil.h"
#include "Walls/CylindricalWall.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"
#include "Walls/InfiniteWallWithHole.h"
#include "Walls/Screw.h"

#include "Boundaries/AngledPeriodicBoundary.h"
#include "Boundaries/ChuteInsertionBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/HopperInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/CircularPeriodicBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/LeesEdwardsBoundary.h"

#include "Species/LinearViscoelasticSpecies.h"

class FullRestartTest : public DPMBase
{
public:
    
    void setupInitialConditions() override
    {
        LinearViscoelasticSpecies* s = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        
        Coil coil(Vec3D(0.1, 0.2, 0.3), 0.4, 0.5, 0.6, 0.7, 0.8);
        coil.setSpecies(speciesHandler.getObject(0));
        
        CylindricalWall cylindricalWall;
        cylindricalWall.setSpecies(speciesHandler.getObject(0));
        cylindricalWall.set(0.1);
        
        AxisymmetricIntersectionOfWalls axisymmetricIntersectionOfWalls;
        axisymmetricIntersectionOfWalls.setSpecies(speciesHandler.getObject(0));
        axisymmetricIntersectionOfWalls.setPosition(Vec3D(0.1, 0.2, 0.3));
        axisymmetricIntersectionOfWalls.setAxis(Vec3D(0.4, 0.5, 0.6));
        axisymmetricIntersectionOfWalls.addObject(Vec3D(0.7, 0.8, 0.9), 0.11 * Vec3D(0.7, 0.8, 0.9));
        
        IntersectionOfWalls intersectionOfWalls;
        intersectionOfWalls.setSpecies(speciesHandler.getObject(0));
        intersectionOfWalls.addObject(Vec3D(0.1, 0.2, 0.3), Vec3D(0.4, 0.5, 0.6));
        
        InfiniteWall infiniteWall;
        infiniteWall.setSpecies(speciesHandler.getObject(0));
        infiniteWall.set(Vec3D(0.1, 0.2, 0.3), 0.4 * Vec3D(0.1, 0.2, 0.3));
        
        InfiniteWallWithHole infiniteWallWithHole;
        infiniteWallWithHole.setSpecies(speciesHandler.getObject(0));
        infiniteWallWithHole.set(Vec3D(0.1, 0.2, 0.3), 0.4, 0.5);
        
        Screw screw(Vec3D(0.1, 0.2, 0.3), 0.4, 0.5, 0.6, 0.7, 0.8);
        screw.setSpecies(speciesHandler.getObject(0));
        
        wallHandler.copyAndAddObject(coil);
        wallHandler.copyAndAddObject(cylindricalWall);
        wallHandler.copyAndAddObject(axisymmetricIntersectionOfWalls);
        wallHandler.copyAndAddObject(intersectionOfWalls);
        wallHandler.copyAndAddObject(infiniteWall);
        wallHandler.copyAndAddObject(infiniteWallWithHole);
        wallHandler.copyAndAddObject(screw);
        
        AngledPeriodicBoundary angledPeriodicBoundary;
        angledPeriodicBoundary.set(Vec3D(0.1, 0.2, 0.3), Vec3D(0.4, 0.5, 0.6), Vec3D(0.7, 0.8, 0.9));
        
        ChuteInsertionBoundary chuteInsertionBoundary;
        chuteInsertionBoundary.set(new SphericalParticle, 13, Vec3D(0.1, 0.2, 0.3), Vec3D(0.4, 0.5, 0.6), 0.33, 0.44, 0.55,
                                   1.1, 1.2);
        
        CubeInsertionBoundary cubeInsertionBoundary;
        cubeInsertionBoundary.set(new SphericalParticle, 13, Vec3D(0.1, 0.2, 0.3), Vec3D(0.4, 0.5, 0.6),
                                  Vec3D(0.11, 0.22, 0.33), Vec3D(0.44, 0.55, 0.66), 0.77, 0.88);
        
        HopperInsertionBoundary hopperInsertionBoundary;
        hopperInsertionBoundary.set(new SphericalParticle, 13, 0.2, 0.3, 0.66, 0.77, 3.1, 0.69, false, 3, 0.11, 0.21, 0.09,
                                    2.31, 0.001, 30);
        
        PeriodicBoundary periodicBoundary;
        periodicBoundary.set(Vec3D(0.1, 0.2, 0.3), 1.1, 2.2, Vec3D(0.4, 0.5, 0.6));
        
        CircularPeriodicBoundary circularPeriodicBoundary(0.3);
        
        DeletionBoundary deletionBoundary;
        deletionBoundary.set(Vec3D(0.1, 0.2, 0.3), 3.14);
        
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.set(
                [&](double time UNUSED)
                {
                    return 0.1;
                },
                [&](double time UNUSED)
                {
                    return 0.2;
                },
                0.3, 0.4, 0.5, 0.6);
        
        boundaryHandler.copyAndAddObject(angledPeriodicBoundary);
        boundaryHandler.copyAndAddObject(chuteInsertionBoundary);
        boundaryHandler.copyAndAddObject(cubeInsertionBoundary);
        boundaryHandler.copyAndAddObject(hopperInsertionBoundary);
        boundaryHandler.copyAndAddObject(periodicBoundary);
        boundaryHandler.copyAndAddObject(circularPeriodicBoundary);
        boundaryHandler.copyAndAddObject(deletionBoundary);
        boundaryHandler.copyAndAddObject(leesEdwardsBoundary);
        
        SphericalParticle baseParticle;
        baseParticle.setSpecies(speciesHandler.getObject(0));
        particleHandler.copyAndAddObject(baseParticle);
    }
};

int main()
{
    FullRestartTest normal;
    normal.setName("FullRestartSelfTest");
    normal.setupInitialConditions();
    normal.writeRestartFile();
    normal.restartFile.setFileType(FileType::MULTIPLE_FILES);
    normal.setName("FullRestartWithCounterSelfTest");
    std::cout << normal.restartFile.getFullName();
    normal.writeRestartFile();
    
    FullRestartTest restart;
    restart.readRestartFile("FullRestartSelfTest.restart");
    restart.setName("FullRestartSelfTestRestarted");
    restart.writeRestartFile();
    restart.write(std::cout);
    
    restart.readRestartFile("FullRestartWithCounterSelfTest.restart.1");
    restart.setName("FullRestartWithCounterSelfTestRestarted");
    restart.writeRestartFile();
}
