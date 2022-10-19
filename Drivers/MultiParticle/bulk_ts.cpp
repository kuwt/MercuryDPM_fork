//Copyright (c) 2013-2022, The MercuryDPM Developers Team. All rights reserved.
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

// Example 1 - Single clump in the box

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/MultiParticle.h"
#include "clump/clump_io.h"
#include "clump/mercury3Dclump.h"
# include <stdlib.h>

Mdouble f_min = -15; Mdouble f_max = 15;
Mdouble f_minz = -100; Mdouble f_maxz = 100;
Mdouble margin = 5;

class multiParticleT1 : public Mercury3Dclump
{
public:
    explicit  multiParticleT1()
    {
        setGravity(Vec3D(0.0, 0.0, -10.0));
        setName("bulk_ts");
        setXBallsAdditionalArguments("-solidf -v0");
        setXMax(f_max);
        setYMax(f_max);
        setZMax(f_maxz);
        setXMin(f_min);
        setYMin(f_min);
        setZMin(f_minz);
        load_clumps(data);
        clump_mass = data.mass[clump_index];
    }

    void setClumpDamping(Mdouble damp){ clump_damping = damp;}

    void setClumpIndex(Mdouble index){ clump_index = index;}

    Mdouble getClumpMass(){return clump_mass;}

    void setupInitialConditions() override
    {
        // Generate random distribution of clumps in the bounding box
        setClumpIndex(1);

        int N_clumps = 100;
        for (int particle = 0; particle < N_clumps; particle ++)
        {
            MultiParticle p0;
            p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to MultiParticle 1
            p0.setMaster();
            clump_data dataR = rotate_clump(data, clump_index, uniform_random_pds()); // Rotate clump arbitrarily
            Vec3D pos = Vec3D(f_min + margin + random_double(f_max - f_min - 2*margin),
                              f_min + margin + random_double(f_max - f_min - 2*margin),
                              f_minz + margin + random_double(f_maxz - f_minz - 2*margin));
            //pos = Vec3D(0.1 * particle,1 * particle,5 * particle);
            p0.setPosition(pos);
            for (int j = 0; j < dataR.pebbles_r[clump_index].size(); j++) {
                p0.addSlave(Vec3D(dataR.pebbles_x[clump_index][j],
                                  dataR.pebbles_y[clump_index][j],
                                  dataR.pebbles_z[clump_index][j]),
                            dataR.pebbles_r[clump_index][j]);
            }

            p0.setPrincipalDirections(
                    Matrix3D(dataR.pd[clump_index][0], dataR.pd[clump_index][1], dataR.pd[clump_index][2],
                             dataR.pd[clump_index][3], dataR.pd[clump_index][4], dataR.pd[clump_index][5],
                             dataR.pd[clump_index][6], dataR.pd[clump_index][7], dataR.pd[clump_index][8]));

            p0.setInitInertia(
                    MatrixSymmetric3D(dataR.toi[clump_index][0], dataR.toi[clump_index][1], dataR.toi[clump_index][2],
                                      dataR.toi[clump_index][4], dataR.toi[clump_index][5],
                                      dataR.toi[clump_index][8]));
            p0.setMassMultiparticle(dataR.mass[clump_index]);
            double maxVel = 1; double minVel = -1;
            double maxAVel = 5; double minAVel = -5;

            Vec3D Vel = Vec3D(minVel + random_double(maxVel-minVel),
                              minVel + random_double(maxVel-minVel),
                              minVel + random_double(maxVel-minVel));

            Vec3D AVel = Vec3D(minAVel + random_double(maxAVel-minAVel),
                               minAVel + random_double(maxAVel-minAVel),
                               minAVel + random_double(maxAVel-minAVel));

            p0.setVelocity(Vel);
            p0.setAngularVelocity(AVel);
            particleHandler.copyAndAddObject(p0);
        }


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
    clump_data data;
    Mdouble clump_mass;
    Mdouble clump_damping = 10;
};

int main(int argc, char* argv[])
{
    multiParticleT1 problem;
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(0.0); // sets the species type-0 density
    species->setDissipation(0.0);
    species->setStiffness(1e6);
    const Mdouble collisionTime = species->getCollisionTime(problem.getClumpMass());
    problem.setClumpDamping(10);
    problem.setTimeStep(collisionTime / 50.0);
    problem.setSaveCount(400);
    problem.setTimeMax(2);
    problem.removeOldFiles();
    problem.solve();
    return 0;
}


