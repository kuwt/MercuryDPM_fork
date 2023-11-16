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

// Example 2 - T-bar featuring Dzhanibekov effect

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/ClumpParticle.h"
#include "../ClumpHeaders/ClumpInput.h"
#include "../ClumpHeaders/Mercury3DClump.h"
# include <stdlib.h>

Mdouble f_min = -10; Mdouble f_max = 10;

class ChangingTOIParticle : public Mercury3Dclump
{
public:
    explicit  ChangingTOIParticle()
    {
        setGravity(Vec3D(0.0, 0.0, 0.0));
        setName("TBar");
        setXBallsAdditionalArguments("-solidf -v0");
        setXMax(f_max);
        setYMax(f_max);
        setZMax(f_max);
        setXMin(f_min);
        setYMin(f_min);
        setZMin(f_min);
        LoadClumps(data);
        setClumpIndex(0);
        clump_mass = data.mass[clump_index];
    }

    void setClumpDamping(Mdouble damp){ clump_damping = damp;}

    void setClumpIndex(Mdouble index){ clump_index = index;}

    Mdouble getClumpMass(){return clump_mass;}

    void setupInitialConditions() override
    {
        // Generate single clump
        setClumpIndex(1);
        ClumpParticle p0;
        p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to Clump 1
        p0.setClump();
        DoubleVector urpds = {0.707, 0.707, 0, -0.707, 0.707, 0, 0, 0, 1};
        // DoubleVector urpds = UniformRandomPDs();
        Vec3D angVel = 2 * Vec3D( urpds[0], urpds[1], urpds[2]);
        data = RotateClump(data, clump_index, urpds); // Rotate clump arbitrarily
        p0.setRadius(data.pebblesR[clump_index][0]);
        Vec3D pos = Vec3D(0, 0, 0);
        p0.setPosition(pos);

        //p0.setInitPrincipalDirections(
        //        Matrix3D(data.pd[clump_index][0], data.pd[clump_index][1], data.pd[clump_index][2],
        //                 data.pd[clump_index][3], data.pd[clump_index][4], data.pd[clump_index][5],
        //                 data.pd[clump_index][6], data.pd[clump_index][7], data.pd[clump_index][8]));

        p0.setPrincipalDirections(
                    Matrix3D(data.pd[clump_index][0], data.pd[clump_index][1], data.pd[clump_index][2],
                             data.pd[clump_index][3], data.pd[clump_index][4], data.pd[clump_index][5],
                             data.pd[clump_index][6], data.pd[clump_index][7], data.pd[clump_index][8]));
        p0.setInitInertia(
                    MatrixSymmetric3D(data.toi[clump_index][0], data.toi[clump_index][1], data.toi[clump_index][2],
                                      data.toi[clump_index][4], data.toi[clump_index][5],
                                      data.toi[clump_index][8]));


        for (int j = 0; j < data.pebblesR[clump_index].size(); j++) {
            p0.addPebble(Vec3D(data.pebblesX[clump_index][j],
                               data.pebblesY[clump_index][j],
                               data.pebblesZ[clump_index][j]),
                         data.pebblesR[clump_index][j]);
        }

        std::cout<<"CLUMP MASS set = "<<data.mass[clump_index]<<std::endl;
        p0.setClumpMass(data.mass[clump_index]);
        p0.setAngularVelocity(angVel);
        
        std::cout<<"CLUMP MASS get = "<<p0.getMass()<<std::endl;
        p0.setDamping(clump_damping);
        particleHandler.copyAndAddObject(p0);


        // Rectangular box
        wallHandler.clear();
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0, 0, getZMax()));
        wallHandler.copyAndAddObject(w0);
    }
private:
    int clump_index;
    ClumpData data;
    Mdouble clump_mass;
    Mdouble clump_damping = 10;
};

int main(int argc, char* argv[])
{
    ChangingTOIParticle problem;
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0); // sets the species type-0 density
    //species->setConstantRestitution(0);
    std::cout<<species->getConstantRestitution()<<std::endl;
    species->setDissipation(0.0);
    species->setStiffness(1e6);
    const Mdouble collisionTime = species->getCollisionTime(problem.getClumpMass());
    problem.setClumpDamping(0);
    problem.setTimeStep(collisionTime / 50.0);

    // Quick demonstration
    problem.setSaveCount(5000);
    problem.setTimeMax(100);


    // For time averaging featuring equipartition - takes a while
    // problem.setSaveCount(50000);
    // problem.setTimeMax(100000);

    problem.removeOldFiles();
    problem.solve();
    return 0;
}
