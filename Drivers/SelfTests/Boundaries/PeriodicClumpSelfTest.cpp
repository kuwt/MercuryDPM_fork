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

// Self-test -- single clump in a periodic box

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/MultiParticle.h"
#include "../../MultiParticle/clump/ClumpIO.h"
#include "../../MultiParticle/clump/Mercury3DClump.h"
# include <stdlib.h>
#include "Boundaries/PeriodicBoundary.h"

Mdouble f_min = -4; Mdouble f_max = 4;

class clumpTest : public Mercury3Dclump
{
public:
    explicit  clumpTest()
    {
        setGravity(Vec3D(0.0, 0.0, 0.0));
        setName("periodicClumpSelfTest");
        setXBallsAdditionalArguments("-solidf -v0");
        setXMax(f_max);
        setYMax(f_max);
        setZMax(f_max);
        setXMin(f_min);
        setYMin(f_min);
        setZMin(f_min);
        load_clumps(data);
        setClumpIndex(0);
        clump_mass = data.mass[clump_index];
    }

    void setClumpDamping(Mdouble damp){ clump_damping = damp;}

    void setClumpIndex(Mdouble index){ clump_index = index;}

    Mdouble getClumpMass(){return clump_mass;}

    void setupInitialConditions() override
    {
    
        // Periodic box
        auto per_x = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_x->set(Vec3D(1, 0, 0), getXMin(), getXMax());

        auto per_y = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_y->set(Vec3D(0, 1, 0), getYMin(), getYMax());

        auto per_z = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_z->set(Vec3D(0, 0, 1), getZMin(), getZMax());
    
    
        // Generate single clump
        setClumpIndex(1);
        MultiParticle p0;
        p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to MultiParticle 1
        p0.setMaster();
        p0.setRadius(data.pebbles_r[clump_index][0]);
        Vec3D pos = Vec3D(0, 0, 0);
        p0.setPosition(pos);
        for (int j = 0; j < data.pebbles_r[clump_index].size(); j++) {
            p0.addSlave(Vec3D(data.pebbles_x[clump_index][j],
                                  data.pebbles_y[clump_index][j],
                                  data.pebbles_z[clump_index][j]),
                            data.pebbles_r[clump_index][j]);
        }
        p0.setPrincipalDirections(
                    Matrix3D(data.pd[clump_index][0], data.pd[clump_index][1], data.pd[clump_index][2],
                             data.pd[clump_index][3], data.pd[clump_index][4], data.pd[clump_index][5],
                             data.pd[clump_index][6], data.pd[clump_index][7], data.pd[clump_index][8]));
        p0.setInitInertia(
                    MatrixSymmetric3D(data.toi[clump_index][0], data.toi[clump_index][1], data.toi[clump_index][2],
                                      data.toi[clump_index][4], data.toi[clump_index][5],
                                      data.toi[clump_index][8]));
        p0.setMassMultiparticle(data.mass[clump_index]);
        p0.setAngularVelocity(Vec3D(0,20,0.01));
        p0.setVelocity(Vec3D(8,0,0));
	p0.setDamping(clump_damping);
        particleHandler.copyAndAddObject(p0);

    }
private:
    int clump_index;
    clump_data data;
    Mdouble clump_mass;
    Mdouble clump_damping = 0;
};

int main(int argc, char* argv[])
{
    clumpTest problem;
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
    problem.setSaveCount(100);
    problem.setTimeMax(1);



    problem.removeOldFiles();
    problem.solve();

    MultiParticle* p = dynamic_cast<MultiParticle*>(problem.particleHandler.getLastObject());
    Vec3D angVel = p->getAngularVelocity();
    Vec3D known_angVel = Vec3D(-1.30004, 19.9932, 0.944699);

    Vec3D pos = p->getPosition();
    Vec3D known_pos = Vec3D(-0.201397, -2.37542, -0.284277);


    helpers::check(angVel,known_angVel,1e-4, "Angular velocity check");
    helpers::check(pos,known_pos,1e-4, "Pebble position");
    
    return 0;
}
