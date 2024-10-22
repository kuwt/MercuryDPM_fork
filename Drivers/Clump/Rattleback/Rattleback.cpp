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

// Example - Rattleback - oblique spinner that switches the direction of rotation when pushed the "wrong" direction

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/ClumpParticle.h"
#include "../ClumpHeaders/ClumpInput.h"
#include "../ClumpHeaders/Mercury3DClump.h"
# include <stdlib.h>

Mdouble f_min = -100; Mdouble f_max = 100;
int SAVECOUNT = 400;

class ChangingTOIParticle : public Mercury3Dclump
{
public:
    explicit  ChangingTOIParticle()
    {
        setGravity(Vec3D(0.0, 0.0, -10.0));
        setName("Rattleback");
        setXBallsAdditionalArguments("-solidf -v0");
        setXMax(2*f_max);
        setYMax(2*f_max);
        setZMax(f_max);
        setXMin(2*f_min);
        setYMin(2*f_min);
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
        // Generate rattleback
        setClumpIndex(3);
        ClumpParticle p0;
        p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to Clump 1
        p0.setClump();
        //DoubleVector upds = UniformRandomPDs();
        DoubleVector upds = {1, 0, 0, 0, 0, 1, 0, -1, 0};
        data = RotateClump(data, clump_index, upds); // here you can try different seeds
        p0.setRadius(data.pebblesR[clump_index][0]);
        Vec3D pos = Vec3D(0, 0, -94);
        p0.setPosition(pos);
        for (int j = 0; j < data.pebblesR[clump_index].size(); j++) {
            p0.addPebble(Vec3D(data.pebblesX[clump_index][j],
                               data.pebblesY[clump_index][j],
                               data.pebblesZ[clump_index][j]),
                         data.pebblesR[clump_index][j]);
        }
        p0.setPrincipalDirections(
                Matrix3D(data.pd[clump_index][0], data.pd[clump_index][1], data.pd[clump_index][2],
                         data.pd[clump_index][3], data.pd[clump_index][4], data.pd[clump_index][5],
                         data.pd[clump_index][6], data.pd[clump_index][7], data.pd[clump_index][8]));
        p0.setInitInertia(
                MatrixSymmetric3D(data.toi[clump_index][0], data.toi[clump_index][1], data.toi[clump_index][2],
                                  data.toi[clump_index][4], data.toi[clump_index][5],
                                  data.toi[clump_index][8]));
        p0.setClumpMass(data.mass[clump_index]);
        double mag = 0;
        p0.setAngularVelocity(Vec3D(0,0,-0.7));
        p0.setVelocity(Vec3D(0,0,0));

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

    void actionsAfterTimeStep() override {
        if (getNumberOfTimeSteps()/SAVECOUNT == getNumberOfTimeSteps()/(double) SAVECOUNT ){
            for (std::vector<BaseParticle*>::iterator it= particleHandler.begin(); it!=particleHandler.end(); ++it){
                if ((*it)->isClump()) {

                    std::cout<<"Saving timestep "<<getNumberOfTimeSteps() <<std::endl;

                    // Add velocity to log
                    std::ofstream funct("ClumpSeq.txt", std::ios_base::app | std::ios_base::out);
                    funct << getNumberOfTimeSteps() << " "
                          << static_cast<ClumpParticle*>(*it)->getPosition() << " "
                          << static_cast<ClumpParticle*>(*it)->getPrincipalDirections_e1() << " "
                          << static_cast<ClumpParticle*>(*it)->getPrincipalDirections_e2() << " "
                          << static_cast<ClumpParticle*>(*it)->getPrincipalDirections_e3() << " "
                    <<"\n";
                    funct.close();

                }
            }

        }
    }


private:
    int clump_index;
    ClumpData data;
    Mdouble clump_mass;
    Mdouble clump_damping = 1;
};

int main(int argc, char* argv[])
{
    ChangingTOIParticle problem;
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0); // sets the species type-0 density
    species->setDissipation(30.0);
    species->setStiffness(1e6);

    species->setSlidingFrictionCoefficient(0.7);
    species->setSlidingStiffness(5e5);
    species->setRollingFrictionCoefficient(0.0);
    species->setRollingStiffness(5e5);

    const Mdouble collisionTime = species->getCollisionTime(problem.getClumpMass());
    problem.setClumpDamping(0.01);
    problem.setTimeStep(collisionTime / 50.0);
    problem.setSaveCount(SAVECOUNT);
    problem.setHGridMaxLevels(1);
    problem.setTimeMax(100.0);
    problem.removeOldFiles();
    problem.solve();
    return 0;
}
