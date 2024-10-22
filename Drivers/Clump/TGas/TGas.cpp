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

// Example - Clumps in the periodic box

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/ClumpParticle.h"
#include "../ClumpHeaders/ClumpInput.h"
#include "../ClumpHeaders/Mercury3DClump.h"
#include <stdlib.h>
#include <CMakeDefinitions.h>

Mdouble f_min = -20; Mdouble f_max = 20; // Size of the box and the margin/clearance for clump seeds
Mdouble margin = 0;


Mdouble av_min = -0; // range of angular velocities
Mdouble av_max = 0;

Mdouble tv_min = -100; // range of translational velocities
Mdouble tv_max = 100;

int N_att = 2000;   // Number of attempts to add particle
std::vector <int> D_h; // log of the number of Dzhanibekov particles

class ChangingTOIParticle : public Mercury3Dclump
{
public:
    explicit  ChangingTOIParticle()
    {
        setGravity(Vec3D(0.0, 0.0, -0.0));
        setName("TGas");
        setXBallsAdditionalArguments("-solidf -v0");
        setXMax(f_max);
        setYMax(f_max);
        setZMax(f_max); // Unbounded domain
        setXMin(f_min);
        setYMin(f_min);
        setZMin(f_min);
        LoadClumps(data);
        clump_mass = data.mass[clump_index];

    }

    void setClumpDamping(Mdouble damp){ clump_damping = damp;}

    void setClumpIndex(Mdouble index){ clump_index = index;}

    Mdouble getClumpMass(){return clump_mass;}

    void setupInitialConditions() override
    {
        setParticlesWriteVTK(1);

        /* Double periodic + bottom wall + unlimited top
        auto per_x = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_x->set(Vec3D(1, 0, 0), getXMin(), getXMax());
        auto per_y = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_y->set(Vec3D(0, 1, 0), getYMin(), getYMax());
        wallHandler.clear();
        InfiniteWall w0;
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        */

        /* Rectangular box
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
        */


        // Periodic box
        auto per_x = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_x->set(Vec3D(1, 0, 0), getXMin(), getXMax());

        auto per_y = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_y->set(Vec3D(0, 1, 0), getYMin(), getYMax());

        auto per_z = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_z->set(Vec3D(0, 0, 1), getZMin(), getZMax());


        /*
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(1); // sets particle radius
        p0.setPosition(Vec3D(0., 0., 0.)); // sets particle position
        p0.setVelocity(Vec3D(0., 0., 0.));// sets particle velocity
        particleHandler.copyAndAddObject(p0);
        */


        // Generate a dense packing of Clumps
        setClumpIndex(1);
        int N_created = 0;
        for (int part = 0; part<N_att; part++) {

                ClumpParticle p0;
                p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to Clump 1
            p0.setClump();
                ClumpData rdata = RotateClump(data, clump_index, UniformRandomPDs()); // Rotate clump arbitrarily



                p0.setRadius(rdata.pebblesR[clump_index][0]);

                for (int j = 0; j < rdata.pebblesR[clump_index].size(); j++) {
                    p0.addPebble(Vec3D(rdata.pebblesX[clump_index][j],
                                       rdata.pebblesY[clump_index][j],
                                       rdata.pebblesZ[clump_index][j]),
                                 rdata.pebblesR[clump_index][j]);
                }
                p0.setPrincipalDirections(
                        Matrix3D(rdata.pd[clump_index][0], rdata.pd[clump_index][1], rdata.pd[clump_index][2],
                                 rdata.pd[clump_index][3], rdata.pd[clump_index][4], rdata.pd[clump_index][5],
                                 rdata.pd[clump_index][6], rdata.pd[clump_index][7], rdata.pd[clump_index][8]));
                p0.setInitInertia(
                        MatrixSymmetric3D(rdata.toi[clump_index][0], rdata.toi[clump_index][1], rdata.toi[clump_index][2],
                                          rdata.toi[clump_index][4], rdata.toi[clump_index][5],
                                          rdata.toi[clump_index][8]));
            p0.setClumpMass(rdata.mass[clump_index]);

                p0.setDamping(clump_damping);


                Vec3D pos = Vec3D(f_min + margin + RandomDouble(f_max - f_min - 2 * margin),
                                  f_min + margin + RandomDouble(f_max - f_min - 2 * margin),
                                  f_min + margin + RandomDouble(f_max - f_min - 2 * margin));

                p0.setPosition(pos);

                Vec3D angVel = Vec3D(av_min + RandomDouble(av_max - av_min),
                                     av_min + RandomDouble(av_max - av_min),
                                     av_min + RandomDouble(av_max - av_min));

                Vec3D vel = Vec3D(tv_min + RandomDouble(tv_max - tv_min),
                                  tv_min + RandomDouble(tv_max - tv_min),
                                  tv_min + RandomDouble(tv_max - tv_min));

                // No motion with zero initial conditions
                //angVel = Vec3D(0,0,0);
                //vel = Vec3D(0,0,0);

                p0.setAngularVelocity(angVel);
                p0.setVelocity(vel);


                if (checkClumpForInteractionPeriodic(p0)) {
                    particleHandler.copyAndAddObject(p0);
                    N_created++;
                }
            }

        std::cout<<"Number of particles created: "<<N_created<<std::endl;

    }

    void actionsAfterTimeStep() override {
        int D_num = 0;
        for (std::vector<BaseParticle*>::iterator it= particleHandler.begin(); it!=particleHandler.end(); ++it){
            if ((*it)->isClump()) {
                D_num += (int) static_cast<ClumpParticle*>(*it)->getDzhanibekovParticle();
            }
        }
        D_h.push_back(D_num);
    }
private:
    int clump_index;
    ClumpData data;
    Mdouble clump_mass;
    Mdouble clump_damping = 0;
};



int main(int argc, char* argv[])
{
    ChangingTOIParticle problem;
    problem.setNumberOfOMPThreads(helpers::readFromCommandLine(argc, argv, "-omp",4));
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0); // sets the species type-0 density
    //species->setConstantRestitution(0);
    std::cout<<species->getConstantRestitution()<<std::endl;
    species->setDissipation(0.0);
    species->setStiffness(1e6);
    const Mdouble collisionTime = species->getCollisionTime(problem.getClumpMass());
    problem.setClumpDamping(0);
    std::cout <<"coll time"<<collisionTime<<std::endl;
    problem.setTimeStep(collisionTime / 50.0);

    // Quick demonstration
    problem.setSaveCount(50);
    problem.setTimeMax(0.1);


    problem.removeOldFiles();
    problem.solve();


    // Return the log of the Dzh particles
    std::ofstream D_hyst;  D_hyst.open ("Dzh_hystory.txt");
    for (int i = 0; i<D_h.size(); i+=50){D_hyst <<D_h[i]<<std::endl;}
    D_hyst.close();

    return 0;
}
