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


#include <iostream>
#include "Species/HertzianViscoelasticFrictionChargedBondedSpecies.h"
#include "DPMBase.h"

/**ChargedBondedParticleUnitTest
 * Two particles are bonded together and relaxed. Useful to test the changed bonded particles
 */
class ChargedBondedParticleUnitTest : public DPMBase
{
public:
    void setupInitialConditions() override {
        setXMax(4);
        setYMax(1);
        setZMax(10);
        setSystemDimensions(3);
        setParticleDimensions(3);
        
        SphericalParticle P0, P1;
        
        P0.setSpecies(speciesHandler.getObject(0));
        P1.setSpecies(speciesHandler.getObject(1));
        
        P0.setPosition(Vec3D(1.9, 0.5, 6.0));
        P1.setPosition(Vec3D(2.1, 0.5, 6.0));
        
        P0.setVelocity(Vec3D(0, 0, 0));
        P1.setVelocity(Vec3D(0, 0, 0));
        
        P0.setRadius(0.25);
        P1.setRadius(0.25);
        
        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);
        
        dynamic_cast<ChargedBondedInteraction*>(interactionHandler.getInteraction(particleHandler.getObject(0),
                                                                                  particleHandler.getObject(1),
                                                                                  0))->bond();
    }
};

int main(int argc, char* argv[])
{
    logger(INFO, "Species is neither Sinter nor LinearPlasticViscoelastic");
    
    //Putting details in more user-friendly terms
    //the maximum force exerted, i.e. when particles are in contact
    double maximumForce = 3;
    //the range of the force.
    //Note that the range is measured from the EDGE of the particle!
    double forceRange = 3;
    //based on user inputs, calculated the necessary adhesion stiffness
    double adStiffness = maximumForce / forceRange;
    //the strength of the force holding particles together
    double bondStrength = 100.0;
    
    ChargedBondedParticleUnitTest dpm;
    
    //setting the species of particles
    auto species = dpm.speciesHandler.copyAndAddObject(HertzianViscoelasticFrictionChargedBondedSpecies());
    
    //setting the material properties of the particles
    species->setDensity(6. / constants::pi);
    species->setEffectiveElasticModulus(1000.0);
    species->setSlidingFrictionCoefficient(0.0);
    species->setAdhesionForceMax(maximumForce);
    species->setAdhesionStiffness(adStiffness);
    species->setBondForceMax(bondStrength);
    species->setBondDissipation(0.2);
    species->setVanDerWaalsForceMax(1);
    species->setVanDerWaalsStiffness(10);
    species->setCharge(-1);
    
    //setting a second species of particles
    auto species2 = dpm.speciesHandler.copyAndAddObject(species);
    species2->setCharge(1);
    
    //Giving a name for the output file
    dpm.setName("ChargedBondedParticleUnitTest");
    //setting the time step of the problem
    dpm.setTimeStep(3e-3);
    dpm.setSaveCount(30);
    //setting gravity to zero to ensure only forces acting are inter-particle forces!
    dpm.setGravity(Vec3D(0, 0, 0));
    //setting the duration of the simulation in "simulation seconds" (determined by
    //the dimensions used in setup)
    dpm.setTimeMax(3.);
    //solving the problem!
    dpm.solve(argc, argv);
    
    helpers::writeToFile("ChargedBondedParticleUnitTest.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'total Energy [J]'\n"
                         "p 'ChargedBondedParticleUnitTest.ene' u 1:($2+$3+$4+$5) w l\n"
    );
    
    Mdouble eneElastic = dpm.getElasticEnergy();
    Mdouble eneKinetic = dpm.getKineticEnergy();
    if (fabs(eneElastic + eneKinetic - (4.03257e-05)) >= 1e-6)
    {
        logger(FATAL, "Particles have the wrong total energy. It is % and should be %", eneElastic+eneKinetic,
               5.33881e-6);
    }
}
