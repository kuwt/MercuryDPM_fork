//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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


#include "Mercury2D.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/SubcriticalMaserBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"


class ShiftingMaserBoundarySelfTest : public Mercury2D
{
public:

    void setupInitialConditions() override
    {
        setName("ShiftingMaserBoundarySelfTest");
        dataFile.setFileType(FileType::MULTIPLE_FILES);

        //set species properties: some standard values
        LinearViscoelasticSpecies species;
        species.setDensity(1);
        species.setCollisionTimeAndRestitutionCoefficient( 2e-2, 0.1, M_PI * 0.08 * 0.08 );
        speciesHandler.copyAndAddObject(species);

        //set time and file properties
        setTimeStep(species.getCollisionTime(1) / 50.0);
        setTimeMax(5.0);
        setSaveCount(0.2 / getTimeStep());

        //set domain size
        setMin({0,0,0});
        setMax({5, 50, 1});

        setGravity({0.3, -1, 0});

        auto wall = new InfiniteWall;
        wall->setSpecies(speciesHandler.getLastObject());
        wall->set(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
        wallHandler.copyAndAddObject(wall);
        wall->set(Vec3D(-1, 0, 0), Vec3D(-0.2, 0, 0));
        wallHandler.copyAndAddObject(wall);
        wall->set(Vec3D(1, 0, 0), Vec3D(1.2, 0, 0));
        wallHandler.copyAndAddObject(wall);

        lifted_ = false;

        //Check if particle is copied correctly when moving
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getLastObject());
        p0.setPosition({0,0,0});
        p0.setVelocity({0,0,0});
        p0.setRadius(0.1);

        auto insb = boundaryHandler.copyAndAddObject(new CubeInsertionBoundary());
        insb->set(&p0, 1, Vec3D(0, 0, 0), Vec3D(1, 3, 0), 
                Vec3D(0, -1, 0), Vec3D(0, -1, 0), 0.07, 0.09);


        write(std::cout, false);
    }

    void actionsAfterTimeStep() override
    {
        if (!lifted_ && particleHandler.getVolume() > 2)
        {
            lifted_ = true;

            std::cerr << "Lifting";
            // wallHandler.removeObject(wallHandler.getLastObject()->getId());
            wallHandler.clear();
            auto wall = new InfiniteWall;
            wall->setSpecies(speciesHandler.getLastObject());
            wall->set(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
            wallHandler.copyAndAddObject(wall);
            wall->set(Vec3D(-1, 0, 0), Vec3D(-0.5, 0, 0));
            wallHandler.copyAndAddObject(wall);

            for (auto p : particleHandler)
                p->addVelocity({1, 0, 0});

            boundaryHandler.clear();
            //set the shifting maser boundary
            auto b0 = boundaryHandler.copyAndAddObject(SubcriticalMaserBoundary());
            b0->set(Vec3D(1.0, 0.0, 0.0), Vec3D(0, -1, 0), 0, 1);
            b0->activateMaser();

            auto db = boundaryHandler.copyAndAddObject(DeletionBoundary());
            db->set(Vec3D(1,0,0), 5);

        }

    }

    bool lifted_;

};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    ShiftingMaserBoundarySelfTest maserSelfTest;
    maserSelfTest.solve();
}

