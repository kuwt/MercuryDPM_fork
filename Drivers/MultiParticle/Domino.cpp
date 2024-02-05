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

// Example 1 - Single clump in the box

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/MultiParticle.h"
#include "clump/ClumpIO.h"
#include "clump/Mercury3DClump.h"
#include <stdlib.h>
#include "Boundaries/PeriodicBoundary.h"

struct dominoes {

    Mdouble margin = 20;  // Distance along y from the edge of the box to the center of the first domino
    int N_dom = 10;        // Number of dominoes
    Mdouble R_peb = 0.75; // Radius of the pebbles forming dominoes
    Mdouble R_cue = 0.65;  // Radius of a cue particle
    Mdouble Vel_cue = 1;  // Velocity of the cue particle
    Mdouble S_peb = 1;    // Spacing of pebbles in domino
    Mdouble S_dom = 5;    // Spacing of dominoes
    int m_peb = 1, n_peb = 2, k_peb = 4; // (m,n,k) are numbers of pebbles in (x,y,z) directions correspondingly
    
    Mdouble x_min = 0;
    Mdouble x_max = 2 * margin + N_dom * S_dom;

    Mdouble y_min = -20;   // Box min in z and z directions
    Mdouble y_max =  20;   // Box max in z and z directions

    Mdouble z_min = -(k_peb-0.5) * S_peb - R_peb; // Box min in z and z directions
    Mdouble z_max =  20; // Box max in z and z directions


    Mdouble mass = 1;
    Mdouble I_xx = mass * (1./12.) * S_peb * S_peb * ( n_peb * n_peb + k_peb * k_peb );
    Mdouble I_yy = mass * (1./12.) * S_peb * S_peb * ( m_peb * m_peb + k_peb * k_peb );
    Mdouble I_zz = mass * (1./12.) * S_peb * S_peb * ( m_peb * m_peb + n_peb * n_peb );

};

dominoes D;


class multiParticleT1 : public Mercury3Dclump
{
public:
    explicit  multiParticleT1()
    {
        setGravity(Vec3D(0.0, 0.0, -10.0));
        setName("Domino");
        setXBallsAdditionalArguments("-solidf -v0");
        setXMax(D.x_max);
        setYMax(D.y_max);
        setZMax(D.z_max);
        setXMin(D.x_min);
        setYMin(D.y_min);
        setZMin(D.z_min);
    }

    void setClumpDamping(Mdouble damp){ clump_damping = damp;}

    void setClumpIndex(Mdouble index){ clump_index = index;}

    Mdouble getClumpMass(){return clump_mass;}

    void setupInitialConditions() override
    {

        // Dominoes
        for (int part = 0; part<D.N_dom; part++) {
            MultiParticle p0;
            p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to MultiParticle 1
            p0.setMaster();
            p0.setRadius(D.R_peb);

            for (int i = 0; i < 2*D.m_peb; i++) {
                for (int j = 0; j < 2*D.n_peb; j++) {
                    for (int k = 0; k < 2*D.k_peb; k++) {

                        p0.addSlave(Vec3D( - D.S_peb * D.m_peb + (i+0.5)*D.S_peb,
                                           - D.S_peb * D.n_peb + (j+0.5)*D.S_peb,
                                           - D.S_peb * D.k_peb + (k+0.5)*D.S_peb), D.R_peb);

            }}}
            p0.setPrincipalDirections(
                    Matrix3D(1, 0, 0,
                             0, 1, 0,
                             0, 0, 1));

            std::cout<<"I_xx = "<<D.I_xx<<std::endl;
            std::cout<<"I_yy = "<<D.I_yy<<std::endl;


            p0.setInitInertia(
                    MatrixSymmetric3D(D.I_xx,       0,      0,
                                                     D.I_yy,      0,
                                                                  D.I_zz));
            p0.setMassMultiparticle(1*D.mass);

            p0.setDamping(clump_damping);



            Vec3D pos = Vec3D(D.margin + part*D.S_dom,
                            0,
                            0);



            p0.setPosition(pos);

            Vec3D angVel(0,0,0);
            Vec3D vel(0,0,0);

            p0.setAngularVelocity(angVel);
            p0.setVelocity(vel);
            particleHandler.copyAndAddObject(p0);
        }
        // Cue
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(D.R_cue); // sets particle radius
        p0.setPosition(Vec3D(0.9*D.margin, 0, D.S_peb*(D.k_peb-0.5))); // sets particle position
        p0.setVelocity(Vec3D(D.Vel_cue, 0., 0.));// sets particle velocity
        particleHandler.copyAndAddObject(p0);

        // Rectangular box
        wallHandler.clear();
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(1));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));
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

    // Domino species
    auto species0 = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());

    // Wall species
    auto species1 = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());

    // Mixed
    auto species01 = problem.speciesHandler.getMixedObject(species0, species1);



    species0->setDensity(1.0); // sets the species type-0 density
    //species->setConstantRestitution(0);
    species0->setSlidingFrictionCoefficient(0.0);
    species0->setSlidingStiffness(5e5);
    species0->setRollingFrictionCoefficient(0.0);
    species0->setRollingStiffness(5e5);
    species0->setDissipation(1.0);
    species0->setStiffness(1e6);
    const Mdouble collisionTime = species0->getCollisionTime(D.mass);

    species01->setSlidingFrictionCoefficient(0.2);
    species01->setSlidingStiffness(5e5);
    species01->setRollingFrictionCoefficient(0.1);
    species01->setDissipation(1.0);
    species01->setRollingStiffness(5e5);
    species01->setStiffness(1e6);



    problem.setClumpDamping(0);
    problem.setTimeStep(collisionTime / 50.0);

    // Quick demonstration
    problem.setSaveCount(500);
    problem.setTimeMax(20);


    // For time averaging featuring equipartition - takes a while
    // problem.setSaveCount(50000);
    // problem.setTimeMax(100000);

    problem.removeOldFiles();
    problem.solve();
    return 0;
}