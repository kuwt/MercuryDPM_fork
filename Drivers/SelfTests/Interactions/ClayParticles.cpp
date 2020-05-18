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
//#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Species/LinearViscoelasticFrictionChargedBondedSpecies.h"
#include "DPMBase.h"
#include "Walls/InfiniteWall.h"
#include "Logger.h"

/// In this file, the rolling behaviour of the tangential spring is tested. This is done by placing one normal partilce on top of a fixed partilce and letting graviry roll it over the other particle until it loses contact.
class ChargedBondedParticleUnitTest : public DPMBase {
public:
    void setupInitialConditions() override {

        //***********************************************Setting System Size and Dimensionality**************************************************
        setXMax(10);
        setYMax(1);
        setZMax(10);
        setXMin(0);
        setYMin(0);
        setZMin(0);
        setSystemDimensions(3);
        setParticleDimensions(3);

        //**************************************Creating Particles and Assigning Species and Properties*******************************************

        //Creating base particles
        SphericalParticle P0,P1,P2,P3,P4;

        //setting the number of conglomerate particles to create
        int nParticles = 6;

        //setting the radius of particles
        double rad = 0.25;

        //Setting the initial position of the first particle placed
        Vec3D basePosition(1.0,0.5,6.0);

        //setting the separation of individual conglomerate particles
        Vec3D macroSeparation(1.25, 0.0, -0.5);

        //setting the initial overlap of bonded particles forming a single
        //conglomerate particle
        auto species = dynamic_cast<LinearViscoelasticFrictionChargedBondedSpecies*>(speciesHandler.getObject(0));
        double delta = species->getBondForceMax()/species->getStiffness();
        Vec3D microSeparationInitial(rad+rad-delta,0.0,0.0);

        //setting the initial velocities of particles
        Vec3D vInitial(2.,0.0,-2.0);

        //creating a counter variable to keep track of the objects being created
        unsigned int ob = 0;

        //a simple example loop to build rod-like particles with positively
        //charged 'ends' and negatively charged 'bodies'
        for (int i = 0; i < nParticles; i++) {
            //assigning species to each member of the conglomerate so that
            //differing charges can be distributed along its length
            P0.setSpecies(speciesHandler.getObject(0));
            P1.setSpecies(speciesHandler.getObject(1));
            P2.setSpecies(speciesHandler.getObject(1));
            P3.setSpecies(speciesHandler.getObject(1));
            P4.setSpecies(speciesHandler.getObject(0));

            Vec3D microSeparation;

            //setting the particles relative to one another in the desired
            //configuration
            //setting initial position in system
            P0.setPosition(Vec3D(basePosition + (i * macroSeparation) + microSeparation));
            microSeparation += microSeparationInitial;
            P1.setPosition(Vec3D(basePosition + (i * macroSeparation) + microSeparation));
            microSeparation += microSeparationInitial;
            P2.setPosition(Vec3D(basePosition + (i * macroSeparation) + microSeparation));
            microSeparation += microSeparationInitial;
            P3.setPosition(Vec3D(basePosition + (i * macroSeparation) + microSeparation));
            microSeparation += microSeparationInitial;
            P4.setPosition(Vec3D(basePosition + (i * macroSeparation) + microSeparation));

            //setting initial velocity, if desired
            P0.setVelocity(Vec3D(vInitial));
            P1.setVelocity(Vec3D(vInitial));
            P2.setVelocity(Vec3D(vInitial));
            P3.setVelocity(Vec3D(vInitial));
            P4.setVelocity(Vec3D(vInitial));

            //setting particle radii
            P0.setRadius(rad);
            P1.setRadius(rad);
            P2.setRadius(rad);
            P3.setRadius(rad);
            P4.setRadius(rad);

            //adding particles to the particle handler
            //(i.e. adding to the system)
            particleHandler.copyAndAddObject(P0);
            particleHandler.copyAndAddObject(P1);
            particleHandler.copyAndAddObject(P2);
            particleHandler.copyAndAddObject(P3);
            particleHandler.copyAndAddObject(P4);



            //Gluing single particles together to make conglomerate particles
            dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(particleHandler.getObject(ob),
                                                                                       particleHandler.getObject(
                                                                                        ob + 1), 0))->bond();
            ob++;
            dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(particleHandler.getObject(ob),
                                                                                       particleHandler.getObject(
                                                                                        ob + 1), 0))->bond();
            ob++;
            dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(particleHandler.getObject(ob),
                                                                                       particleHandler.getObject(
                                                                                        ob + 1), 0))->bond();
            ob++;
            dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(particleHandler.getObject(ob),
                                                                                       particleHandler.getObject(
                                                                                        ob + 1), 0))->bond();
            ob+=2;
        }

        //***************************************Creating Walls and Assigning Species and Properties**********************************************
        InfiniteWall baseWall;
        //setting wall species
        baseWall.setSpecies(speciesHandler.getObject(2));
        //setting the position of the wall
        baseWall.set(Vec3D( 0.0, 0.0,-1.0),Vec3D(0,0,getZMin()));
        wallHandler.copyAndAddObject(baseWall);

        InfiniteWall topWall;
        //setting wall species
        topWall.setSpecies(speciesHandler.getObject(2));
        //setting the position of the wall
        topWall.set(Vec3D( 0.0, 0.0,1.0),Vec3D(0,0,getZMax()));
        wallHandler.copyAndAddObject(topWall);

        InfiniteWall leftWall;
        //setting wall species
        leftWall.setSpecies(speciesHandler.getObject(2));
        //setting the position of the wall
        leftWall.set(Vec3D( -1.0, 0.0,0.0),Vec3D(getXMin(),0,0));
        wallHandler.copyAndAddObject(leftWall);

        InfiniteWall rightWall;
        //setting wall species
        rightWall.setSpecies(speciesHandler.getObject(2));
        //setting the position of the wall
        rightWall.set(Vec3D( 1.0, 0.0,0.0),Vec3D(getXMax(),0,0));
        wallHandler.copyAndAddObject(rightWall);

        InfiniteWall frontWall;
        //setting wall species
        frontWall.setSpecies(speciesHandler.getObject(2));
        //setting the position of the wall
        frontWall.set(Vec3D( 0.0, -1.0,0.0),Vec3D(0,getYMin(),0));
        wallHandler.copyAndAddObject(frontWall);

        InfiniteWall backWall;
        //setting wall species
        backWall.setSpecies(speciesHandler.getObject(2));
        //setting the position of the wall
        backWall.set(Vec3D( 0.0, 1.0,0.0),Vec3D(0,getYMax(),0));
        wallHandler.copyAndAddObject(backWall);
    }
};

Logger<Log::ERROR> testLogger("chargedParticleForceUnitTest");

int main(int argc, char *argv[])
{
    //********************************************************Setting Interaction Details**********************************************************
    //Putting details in more user-friendly terms
    //the maximum force exerted, i.e. when particles are in contact
    double maximumForce = 3;
    //the range of the force.
    //Note that the range is measured from the EDGE of the particle!
    double forceRange = 3;
    //based on user inputs, calculated the necessary adhesion stiffness
    double adStiffness = maximumForce / forceRange;
    //based on user inputs, calculated the necessary adhesion stiffness
    double stiffness = 1600 * adStiffness;
    //the strength of the force holding particles together
    //(Affects separation of bonded particles)
    double bondStrength = 333.0;


    ChargedBondedParticleUnitTest chargedParticleForceUnitTestProblem;
    //setting the species of particles
    auto species1 = chargedParticleForceUnitTestProblem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionChargedBondedSpecies());
    //adding second species to allow for attractive interactions
    auto species2 = chargedParticleForceUnitTestProblem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionChargedBondedSpecies());

    //setting species for system walls
    auto speciesW = chargedParticleForceUnitTestProblem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionChargedBondedSpecies());

    //*************************************ASSIGNING PROPERTIES TO PARTICLE SPECIES*****************************************
    //SPECIES 1
    //setting the material properties of the particles
    species1->setDensity(6./constants::pi);
    species1->setStiffness(stiffness);
    species1->setSlidingFrictionCoefficient(1);
    species1->setSlidingStiffness(2.0/7.0*stiffness);
    species1->setRollingFrictionCoefficient(1);
    species1->setRollingStiffness(2.0/5.0*stiffness);
    //setting the interaction force properties
    //setting adhesion stiffness and force max as required to give the
    //user-requested force range.
    species1->setAdhesionForceMax(maximumForce);
    species1->setAdhesionStiffness(adStiffness);
    //Setting the charge of a particle:
    //Positive (1) negative (-1) or neutral (0)
    species1->setCharge(1);
    //Assigning bond properties:
    //The maximum force exerted by the 'bond'
    //(determines the size of the overlap)
    species1->setBondForceMax(bondStrength);
    //the dissipation of the interaction -
    //used to damp out any unphysical oscillations in the contacting pair
    species1->setBondDissipation(0.2);

    //SPECIES 2
    //setting the material properties of the particles
    species2->setDensity(6./constants::pi);
    species2->setStiffness(stiffness);
    species2->setSlidingFrictionCoefficient(1);
    species2->setSlidingStiffness(2.0/7.0*stiffness);
    species2->setRollingFrictionCoefficient(1);
    species2->setRollingStiffness(2.0/5.0*stiffness);
    //setting the interaction force properties
    //setting adhesion stiffness and force max as required to give the
    //user-requested force range
    species2->setAdhesionForceMax(maximumForce);
    species2->setAdhesionStiffness(adStiffness);
    //Setting the charge of a particle:
    //Positive (1) negative (-1) or neutral (0)
    species2->setCharge(-1);
    //Assigning bond properties:
    //The maximum force exerted by the 'bond'
    //(determines the size of the overlap)
    species2->setBondForceMax(bondStrength);
    //the dissipation of the interaction -
    //used to damp out any unphysical oscillations in the contacting pair
    species2->setBondDissipation(0.2);


    //FOR MIXED (SPECIES1-SPECIES2) INTERACTIONS
    //adding also the ability to alter the mixed-particle-interaction properties
    auto species1_2 = dynamic_cast<LinearViscoelasticFrictionChargedBondedMixedSpecies*>(chargedParticleForceUnitTestProblem.speciesHandler.getMixedObject(species1->getIndex(),species2->getIndex()));
    //for now, simply setting mixed object parameters equal to those of species 1
    //setting the material properties to the average values of the 2 particle species undergoing interaction
    species1_2->mixAll(species1,species2);
    species1_2->setSlidingFrictionCoefficient(1);
    species1_2->setSlidingStiffness(2.0/7.0*stiffness);
    species1_2->setRollingFrictionCoefficient(1);
    species1_2->setRollingStiffness(2.0/5.0*stiffness);

    //*************************************ASSIGNING PROPERTIES TO WALL SPECIES*********************************************
    //SPECIES W
    //setting the material properties of the particles
    speciesW->setDensity(6./constants::pi);
    speciesW->setStiffness(stiffness);
    speciesW->setSlidingFrictionCoefficient(1);
    speciesW->setSlidingStiffness(2.0/7.0*stiffness);
    speciesW->setRollingFrictionCoefficient(1);
    speciesW->setRollingStiffness(2.0/5.0*stiffness);
    //setting the interaction force properties
    //setting adhesion stiffness and force max as required to give the
    //user-requested force range
    speciesW->setAdhesionForceMax(maximumForce);
    speciesW->setAdhesionStiffness(adStiffness);
    //Setting the charge of a particle:
    //Positive (1) negative (-1) or neutral (0)
    //charge by default set to zero. However, can be changed, e.g. to induce
    //a field throughout the system and encourage alignment
    speciesW->setCharge(0);
    //Assigning bond properties:
    //The maximum force exerted by the 'bond'
    //(determines the size of the overlap)
    //For walls, setting to zero as no need for bonded interactions here.
    speciesW->setBondForceMax(0);
    //the dissipation of the interaction -
    //used to damp out any unphysical oscillations in the contacting pair
    speciesW->setBondDissipation(0.2);


    //FOR MIXED (PARTICLE-WALL) INTERACTIONS
    //walls and species 1 particles
    auto speciesW_1 = dynamic_cast<LinearViscoelasticFrictionChargedBondedMixedSpecies*>(chargedParticleForceUnitTestProblem.speciesHandler.getMixedObject(speciesW->getIndex(),species1->getIndex()));
    //setting the material properties to the average values of the 2 particle species undergoing interaction
    speciesW_1->mixAll(speciesW,species1);
    speciesW_1->setSlidingFrictionCoefficient(1);
    speciesW_1->setSlidingStiffness(2.0/7.0*stiffness);
    speciesW_1->setRollingFrictionCoefficient(1);
    speciesW_1->setRollingStiffness(2.0/5.0*stiffness);

    //walls and species 1 particles
    auto speciesW_2 = dynamic_cast<LinearViscoelasticFrictionChargedBondedMixedSpecies*>(chargedParticleForceUnitTestProblem.speciesHandler.getMixedObject(speciesW->getIndex(),species2->getIndex()));
    //setting the material properties to the average values of the 2 particle species undergoing interaction
    speciesW_2->mixAll(speciesW,species2);
    speciesW_2->setSlidingFrictionCoefficient(1);
    speciesW_2->setSlidingStiffness(2.0/7.0*stiffness);
    speciesW_2->setRollingFrictionCoefficient(1);
    speciesW_2->setRollingStiffness(2.0/5.0*stiffness);

    //*************************************ASSIGNING GENERAL SIMULATION PARAMETERS******************************************
    //Giving a name for the output file
    chargedParticleForceUnitTestProblem.setName("ClayParticles");
    //setting the time step of the problem
    double radius = 0.25;
    double mass = species1->getMassFromRadius(radius);
    double tc = species1->getCollisionTime(mass);
    chargedParticleForceUnitTestProblem.setTimeStep(0.02*tc);
    //setting gravity to zero to ensure only forces acting are inter-particle forces!
    chargedParticleForceUnitTestProblem.setGravity(Vec3D(0,0,-0));
    //setting the duration of the simulation in "simulation seconds" (determined by
    //the dimensions used in setup)
    chargedParticleForceUnitTestProblem.setTimeMax(30.);
    //setting parameters for visualisation using the xballs software package
    chargedParticleForceUnitTestProblem.setXBallsAdditionalArguments("-cmode 8 -solidf -v0");
    //solving the problem!
    chargedParticleForceUnitTestProblem.solve(argc,argv);

    std::vector<BaseParticle*>::iterator pIt = chargedParticleForceUnitTestProblem.particleHandler.begin();

}

























