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

#include <Particles/SphericalParticle.h>
#include "Mercury3D.h"
//#include "Math/ExtendedMath.h"
#include "Species/LinearViscoelasticFrictionReversibleAdhesiveSpecies.h"
//#include "Walls/InfiniteWall.h"
//#include "Walls/IntersectionOfWalls.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"

class MovingWalls : public DPMBase
{
public:
    explicit MovingWalls (std::string name) {
        setName("MovingWalls" + name);
        solve();
        logger(INFO,"%", *wallHandler.getLastObject());
        logger(INFO,"Test successful");
    }

    void setupInitialConditions() override
    {
        //stiffness and distance are set such that all forces are 1
        LinearViscoelasticFrictionReversibleAdhesiveSpecies species;
        species.setStiffnessAndRestitutionCoefficient(100.0,0.5,1.0);
        species.setAdhesionStiffness(100);
        species.setSlidingDissipation(2./7*species.getDissipation());
        species.setSlidingStiffness(2./7*species.getStiffness());
        species.setSlidingFrictionCoefficient(0.5);
        species.setRollingDissipation(2./7*species.getDissipation());
        species.setRollingStiffness(2./7*species.getStiffness());
        species.setRollingFrictionCoefficient(0.5);
        species.setTorsionDissipation(2./7*species.getDissipation());
        species.setTorsionStiffness(2./7*species.getStiffness());
        species.setTorsionFrictionCoefficient(0.5);
        species.setAdhesionStiffness(100);
        species.setAdhesionForceMax(1);
        species.setDensity(6.0 / constants::pi);
        speciesHandler.copyAndAddObject(species);

        Vec3D position = Vec3D(0,0,1);

        //create periodic boundaries around the domain
        if (getName() == "MovingWallsInfiniteWall") {
            InfiniteWall w;
            w.setSpecies(speciesHandler.getObject(0));
            w.set(Vec3D(0,0,1),position);

            w.setAngularVelocity(Vec3D(0,.1,0));
            wallHandler.copyAndAddObject(w);
        } else if (getName() == "MovingWallsIntersectionOfWalls") {
            IntersectionOfWalls w;
            w.setSpecies(speciesHandler.getObject(0));
            w.setPosition(position);
            w.addObject(Vec3D(0,0,1),Vec3D(0,0,0));

            w.setAngularVelocity(Vec3D(0,.1,0));
            wallHandler.copyAndAddObject(w);
        } else if (getName() == "MovingWallsAxisymmetricIntersectionOfWalls") {
            AxisymmetricIntersectionOfWalls w;
            w.setSpecies(speciesHandler.getObject(0));
            w.setPosition(position);
            w.setOrientationViaNormal(Vec3D(0,0,1));
            w.addObject(Vec3D(0,0,1),Vec3D(0,0,0));

            w.setAngularVelocity(Vec3D(0,.1,0));
            wallHandler.copyAndAddObject(w);
        } else {
            logger(ERROR,"WallType not recognized");
        }

        //put particles into the domain in a 3d grid of mesh size distance
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(0.5);
        p.setPosition(position+Vec3D(1,0,-0.49));
        particleHandler.copyAndAddObject(p);

        // define the domain
        setMin(position-Vec3D(1,1,1));
        setMax(position+Vec3D(1,1,1));

        //dataFile.setFileType(FileType::NO_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        setWallsWriteVTK(FileType::MULTIPLE_FILES);
        interactionHandler.setWriteVTK(FileType::MULTIPLE_FILES);
        setParticlesWriteVTK(true);
        setSaveCount(100);
        setTimeStep(0.02*species.getCollisionTime(1.0));
        setTimeMax(20.0*constants::pi);
    }
};



int main()
{
    logger(INFO,"Create wall and attach one particle.");
    MovingWalls problem1("IntersectionOfWalls");
    MovingWalls problem0("InfiniteWall");
    MovingWalls problem2("AxisymmetricIntersectionOfWalls");
    return 0;
}
