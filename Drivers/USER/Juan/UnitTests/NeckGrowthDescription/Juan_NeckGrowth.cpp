//Copyright (c) 2013-2014, The MercuryDPM Developers Team. All rights reserved.
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

#include "Mercury3D.h"
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>

//This UnitTest code reproduces the neck growth between two particles.
class NeckGrowth : public DPMBase{
public:

    void setupInitialConditions() override {

        //++++++++++++++++++++++++++
        //Global parameters:
        setGravity(Vec3D(0.0,0.0,0.0));
        //++++++++++++++++++++++++++

        //++++++++++++++++++++++++++
        //Boundary:
        setXMax(2.0*radius);
        setYMax(radius);
        setZMax(radius);
        setXMin(-getXMax());
        setYMin(-getYMax());
        setZMin(-getZMax());
        setGravity(Vec3D(0.0,0.0,0.0));
        //++++++++++++++++++++++++++

        //++++++++++++++++++++++++++
        //Species of the particles:
        ThermalParticle P0,P1;
        P0.setSpecies(speciesHandler.getObject(0));
        P1.setSpecies(speciesHandler.getObject(0));
        //++++++++++++++++++++++++++

        //++++++++++++++++++++++++++
        //Position and velocity of the particles:
        P0.setPosition({-(1-1e-15)*radius,0,0});
        P0.setRadius(radius);

        P1.setPosition(-P0.getPosition());
        P1.setRadius(radius);

        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);
        //++++++++++++++++++++++++++
    }

    Mdouble radius = 0.5;
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* species = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    NeckGrowth sf;

    //----------------------------------------------
    //Parameters:
    sf.species->setDensity(1005);
    sf.radius = 5e-7;

    const Mdouble restitutionCoefficient = 0.1;
    const Mdouble stiffness = 1e-2 * sf.radius;
    const Mdouble mass = sf.species->getMassFromRadius(sf.radius);

    //----------------------------------------------
    //Compute Dissipation:
    sf.species->setStiffnessAndRestitutionCoefficient(stiffness,restitutionCoefficient,mass);
    //Set elasto-plastic parameters. [k1,k2,kc,phi]:
    sf.species->setPlasticParameters(stiffness, 10*stiffness, stiffness, 0.16);
    //get collision time:
    const Mdouble collisionTime = sf.species->getCollisionTime(mass);
    //----------------------------------------------

    //----------------------------------------------
    sf.species->setSinterType(SINTER_APPROACH::FRENKEL);
//    sf.species->setSinterRate(4e-10/sf.radius);
    sf.species->setSinterAdhesion(0.013*stiffness);
    //adhesiveForce = sinterAdhesion*radius;

    //----------------------------------------------
    //Simulation set up:
    sf.setTimeStep(0.02*collisionTime);
    sf.setTimeMax(150);
    sf.setSaveCount(20000);
    sf.setFileType(FileType::ONE_FILE);
    //----------------------------------------------

    //----------------------------------------------
    //Output file:
    logger(INFO, "Running for Particle radious=%",sf.radius);
    std::string r = helpers::to_string(sf.radius);
    sf.setName("Juan_NeckGrowth"+r);
    sf.writeRestartFile();
    sf.solve();

    //----------------------------------------------
    //Helper [1]:
    std::cout << "Execute 'gnuplot Juan_NeckGrowth+++.gnu' to view output" << std::endl;
    helpers::writeToFile("Juan_NeckGrowth+++.gnu",
                         "set xlabel 'displacement'\n"
                         "set ylabel 'force'\n"
                         "plot 'Juan_NeckGrowth+++.fstat' u 7:9 w lp\n"
    );
    //Helper [2]:
    helpers::writeToFile("Juan_NeckGrowth+++.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'x/a'\n"
                         "plot [0:200] 'Juan_NeckGrowth5e-07.fstat' u ($1):(sqrt($7/5e-07)) w lp lt rgb 'royalblue'");
    //----------------------------------------------
}
