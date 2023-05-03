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
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Logger.h"

/// In this file, the rolling behaviour of the tangential spring is tested. This is done by placing one normal partilce on top of a fixed partilce and letting graviry roll it over the other particle until it loses contact.
class ChargedBondedParticleUnitTest : public DPMBase {
public:
    void setupInitialConditions() override
    {

        //***********************************************Setting System Size and Dimensionality**************************************************
        setXMax(xSize);
        setYMax(ySize);
        setZMax(zSize);
        setXMin(0);
        setYMin(0);
        setZMin(0);
        setSystemDimensions(3);
        setParticleDimensions(3);

        //**************************************Creating Particles and Assigning Species and Properties*******************************************

        //Creating base particles
        SphericalParticle P0,P1,P2,P3,P4;

        //setting the initial overlap of bonded particles forming a single
        //conglomerate particle
        auto species = dynamic_cast<LinearViscoelasticFrictionChargedBondedSpecies*>(speciesHandler.getObject(0));
        double delta = species->getBondForceMax()/species->getStiffness();
        Vec3D microSeparationInitial(radius+radius-delta,0.0,0.0);

        Vec3D microSeparation(0.0,0.0,0.0);

        //setting the overlap between rows of particles forming a hexagon
        double microOverlap = 0.02;

        //setting the initial velocities of particles
        Vec3D vInitial(2.,3.0,-2.0);

        //creating a counter variable to keep track of the objects being created
        int ob = 0;

        //************************************For a simple 2 particle collision*****************************************************************
/*
        P0.setSpecies(speciesHandler.getObject(0));
        P1.setSpecies(speciesHandler.getObject(1));

        P0.setPosition(Vec3D(1.0,0.5,0.5));
        P0.setVelocity(Vec3D(2.0,0,0));
        P0.setRadius(radius);
        particleHandler.copyAndAddObject(P0);


        P1.setPosition(Vec3D(9.0,0.5,0.5));
        P1.setVelocity(Vec3D(-2.0,0,0));
        P1.setRadius(radius);
        particleHandler.copyAndAddObject(P1);

*/

        //****************************************For Hexagons*************************************************************************************
//building a conglomerate particle in a hexagonal shape

        Vec3D basePositionInitial = basePosition;

        for (int n = 0; n < nParticles; n++) {

            //A variable to store the current width of each row
            int nWidthI = nWidthInitial;

            //looping over (half + 1) of the required rows to form the first half of an individual hexagon
            for (int i = 0; i < nRow / 2 + 1; i++) {
                //building an individual row of the hexagon
                for (int j = 0; j < nWidthI; j++) {
                    //choosing species dependent on position within hexagon
                    //putting species 1 on the 'outside' of the particle
                    //checking first if particles are in the first row, in which case they are
                    //all 'external'
                    if (i == 0) {
                        P0.setSpecies(speciesHandler.getObject(1));
                    }
                        //checking next if they are on the far left of the current row
                    else if (j == 0) {
                        P0.setSpecies(speciesHandler.getObject(1));
                    }
                        //or the far right
                    else if (j == (nWidthI - 1)) {
                        P0.setSpecies(speciesHandler.getObject(1));
                    }
                        //and otherwise, setting to species 0
                    else {
                        P0.setSpecies(speciesHandler.getObject(0));
                    }
                    //position, size and velocity are then set irrespective of particle species
                    P0.setPosition(Vec3D(basePosition + microSeparation));
                    P0.setVelocity(Vec3D(vInitial));
                    P0.setRadius(radius);
                    particleHandler.copyAndAddObject(P0);
                    microSeparation += microSeparationInitial;
                }
                //setting up for the next row
                //shifting the base poisition of the next row of particles left by half
                //a particle diameter and down by sqrt(3)
                basePosition -= Vec3D(radius, 0.0, radius * sqrt(3) - microOverlap);
                //increasing the number of particles to be added to the next row
                nWidthI++;
                //resetting the microseparation to zero
                microSeparation.setZero();
            }
            //lowering nWidth in preparation for the next loop...
            nWidthI--;
            //...and, for similar reasons, adjusting the base position
            basePosition += Vec3D(radius, 0.0, radius * sqrt(3) - microOverlap);



            //looping over the remaining rows to form the second half
            //looping over (half + 1) of the required rows to form the first half of an individual hexagon
            for (int i = 0; i < nRow / 2; i++) {
                //decrementing the number of particles in the row as we are now past the widest point
                //of the hexagon!
                nWidthI--;
                //setting up for the next row
                //shifting the base poisition of the next row of particles right by half
                //a particle diameter and down by sqrt(3)
                basePosition += Vec3D(radius, 0.0, -radius * sqrt(3) + microOverlap);
                //resetting the microseparation to zero
                microSeparation.setZero();
                //building an individual row of the hexagon
                for (int j = 0; j < nWidthI; j++) {
                    //choosing species dependent on position within hexagon
                    //putting species 1 on the 'outside' of the particle
                    //checking first if particles are in the last row, in which case they are
                    //all 'external'
                    if (i == (nWidthI - 1)) {
                        P0.setSpecies(speciesHandler.getObject(1));
                    }
                        //checking next if they are on the far left of the current row
                    else if (j == 0) {
                        P0.setSpecies(speciesHandler.getObject(1));
                    }
                        //or the far right
                    else if (j == (nWidthI - 1)) {
                        P0.setSpecies(speciesHandler.getObject(1));
                    }
                        //and otherwise, setting to species 0
                    else {
                        P0.setSpecies(speciesHandler.getObject(0));
                    }


                    P0.setPosition(Vec3D(basePosition + microSeparation));
                    P0.setVelocity(Vec3D(vInitial));
                    P0.setRadius(radius);
                    particleHandler.copyAndAddObject(P0);
                    microSeparation += microSeparationInitial;
                }
            }



            //Gluing single particles together to make conglomerate particles
            //again, for clarity, doing separately for two halves of hexagon

            //resetting the row number
            nWidthI = nWidthInitial;

            //looping over (half -1) of the rows forming the first half of an individual hexagon
            for (int i = 0; i < nRow / 2; i++) {
                //looping again to glue each individual row of the hexagon
                for (int j = 0; j < nWidthI; j++) {
                    //dependent on the position of the particle, there are two possible
                    //manners in which it may need to be glued!
                    //specifically, we treat particles on the right-hand edge as different from others
                    //as they do not have a neighbour on their right
                    //Checking if particle is furthest-right
                    if (j == nWidthI - 1) {
                        //if so, joining only to particles below
                        dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob), particleHandler.getObject(ob + nWidthI), 0))->bond();
                        dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob), particleHandler.getObject(ob + nWidthI + 1), 0))->bond();
                    }
                        //Otherwise, joining also to adjacent particle
                    else {
                        dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob), particleHandler.getObject(ob + nWidthI), 0))->bond();
                        dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob), particleHandler.getObject(ob + nWidthI + 1), 0))->bond();
                        dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob), particleHandler.getObject(ob + 1), 0))->bond();
                    }
                    //moving on to the next particle (object) in the particle handler
                    ob++;
                }
                //Incrementing N width to ensure counting remains correct!
                nWidthI += 1;


            }

            //looping over the remaining rows to complete the gluing process. Gluing this time to the row *above*
            for (int i = 0; i < nRow / 2 + 1; i++) {
                //looping again to glue each individual row of the hexagon
                for (int j = 0; j < nWidthI; j++) {
                    //dependent on the position of the particle, there are two possible
                    //manners in which it may need to be glued!
                    //specifically, we treat particles on the left-hand each edge as different from others
                    //as they do not have two upper neighbours
                    //checking if particle is furthest-left
                    if (j == 0) {
                        //if so, join to particle on the top-right...
                        dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob - nWidthI), particleHandler.getObject(ob), 0))->bond();
                        //and particle to top-left
                        dynamic_cast<ChargedBondedInteraction*>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob - (nWidthI + 1)), particleHandler.getObject(ob),
                                0))->bond();
                    }
                        //Otherwise, joining also to adjacent particle
                    else {
                        //if so, join to particle on the top-right...
                        dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob - nWidthI), particleHandler.getObject(ob), 0))->bond();
                        //and particle to top-left....
                        dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob - (nWidthI + 1)), particleHandler.getObject(ob),
                                0))->bond();
                        //...AND adjacent particle to the left
                        dynamic_cast<ChargedBondedInteraction *>(interactionHandler.getInteraction(
                                particleHandler.getObject(ob - 1), particleHandler.getObject(ob), 0))->bond();

                        std::cout << "Gluing objects " << ob - 1 << " and " << ob << std::endl;
                    }
                    //moving on to the next particle (object) in the particle handler
                    ob++;
                }
                //Incrementing N width to ensure counting remains correct!
                nWidthI -= 1;

            }
            //incrementing the base position to make sure no overlap between particles!
            basePosition = basePositionInitial;
            basePosition += macroSeparation * (n+1);
            microSeparation.setZero();
        }








        //****************************************For simple rods**********************************************************************************
/*
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
            P0.setRadius(radius);
            P1.setRadius(radius);
            P2.setRadius(radius);
            P3.setRadius(radius);
            P4.setRadius(radius);

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
*/
        //***************************************Creating Walls and Assigning Species and Properties**********************************************

        //*******************************************************For Periodic Walls***************************************************************

        if (!solidWalls) {
            PeriodicBoundary zBounds;
            zBounds.set(Vec3D(0.0, 0.0, 1.0), getZMin(), getZMax());
            boundaryHandler.copyAndAddObject(zBounds);

            //setting up static periodic boundaries in the x- and y-directions
            PeriodicBoundary xBounds;
            xBounds.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
            boundaryHandler.copyAndAddObject(xBounds);

            PeriodicBoundary yBounds;
            yBounds.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(yBounds);
        }
        //*******************************************************For Solid Walls******************************************************************
        if (solidWalls) {
            InfiniteWall baseWall;
            //setting wall species
            baseWall.setSpecies(speciesHandler.getObject(2));
            //setting the position of the wall
            baseWall.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));
            wallHandler.copyAndAddObject(baseWall);

            InfiniteWall topWall;
            //setting wall species
            topWall.setSpecies(speciesHandler.getObject(2));
            //setting the position of the wall
            topWall.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0, 0, getZMax()));
            wallHandler.copyAndAddObject(topWall);

            InfiniteWall leftWall;
            //setting wall species
            leftWall.setSpecies(speciesHandler.getObject(2));
            //setting the position of the wall
            leftWall.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0, 0));
            wallHandler.copyAndAddObject(leftWall);

            InfiniteWall rightWall;
            //setting wall species
            rightWall.setSpecies(speciesHandler.getObject(2));
            //setting the position of the wall
            rightWall.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0, 0));
            wallHandler.copyAndAddObject(rightWall);

            InfiniteWall frontWall;
            //setting wall species
            frontWall.setSpecies(speciesHandler.getObject(2));
            //setting the position of the wall
            frontWall.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0, getYMin(), 0));
            wallHandler.copyAndAddObject(frontWall);

            InfiniteWall backWall;
            //setting wall species
            backWall.setSpecies(speciesHandler.getObject(2));
            //setting the position of the wall
            backWall.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0, getYMax(), 0));
            wallHandler.copyAndAddObject(backWall);
        }
    }

    //set and get functions
    //a function to choose if walls bounding the system are solid or periodic
    void setSolidWalls(bool s) {
        solidWalls = s;
    }
    //a function to set the (maximal) strength and range of charge-induced forces
    void setChargeForceAndRange(double F,double R) {
        maxChargeForce = F;
        chargeForceRange = R;
        chargeForceStiffness = F / R;
    }

    //a series of parameters to retrieve the variables used above relating to pseudo-EM forces
    double getMaxChargeForce() {
        return maxChargeForce;
    }
    double getChargeForceRange() {
        return chargeForceRange;
    }
    double getChargeForceStiffness() {
        return chargeForceStiffness;
    }
    //defining a function with which to set the system size
    void setSystemSize(double x, double y, double z) {
        xSize = x;
        ySize = y;
        zSize = z;
    }
    //defining a pair of functions to set and retrieve the radius of an individual spherical 'sub-particle'
    void setSubParticleRadius(double rad) {
        radius = rad;
    }
    double getSubParticleRadius() {
        return radius;
    }

    //a function to set the total number of (conglomerate) particles within the system
    void setParticleNumber(int n) {
        nParticles = n;
    }

    //setting the initial position of the first particle...
    void setBasePosition(Vec3D pos) {
        basePosition = pos;
    }
    //...and a parameter to determine the positions of subsequent particles
    void setMacroSeparation(Vec3D sep) {
        macroSeparation = sep;
    }
    //functions to determine the size and shape of the hexagonal particles formed
    void setRowNumber (int n) {
        nRow = n;
    }
    void setRowSize (int n) {
        nWidthInitial = n;
    }

    //for safety, defining the main variables which affect system behaviour as private
private:
    //a boolean which allows the user to easily choose whether to use periodic boundaries or solid walls
    bool solidWalls;

    //a series of parameters which determine the strength and range of the electromagnetic-style forces which
    //act between particles possessing 'charge'
    //the maximum force due to charge - i.e. that experienced at contact
    double maxChargeForce;
    //the maximal range of the charge-induced force
    double chargeForceRange;
    //the 'stiffness' which determined the two parameters above during calculations
    double chargeForceStiffness;
    //defining a set of parameters with which to set the system size
    double xSize;
    double ySize;
    double zSize;

    //The radius of an individual particle which forms part of a larger, conglomerate particle
    double radius;

    //The number of **full conglomerate particles** to be created in the system
    int nParticles;

    //the position of the first particle put into the system
    Vec3D basePosition;

    //a parameter determining how particles are spread throughout the system
    Vec3D macroSeparation;

    //the size of the 'first row' of a hexagon being created
    int nWidthInitial;
    //The total number of rows forming the hexagon
    int nRow;


};

Logger<Log::ERROR> testLogger("chargedParticleForceUnitTest");

int main(int argc, char *argv[])
{
    //*********************************************************************************************************************************************
    //******************************************************Setting System Details*****************************************************************
    //*********************************************************************************************************************************************

    //******************************************************i) System walls and geometry***********************************************************
    //Creating an instance of the relevant class for this problem
    ChargedBondedParticleUnitTest chargedParticleForceUnitTestProblem;
    //Choosing either a solid-walled system or a system with periodic boundaries
    //true gives solid walls, false gives periodic
    chargedParticleForceUnitTestProblem.setSolidWalls(true);
    //setting the size of the system (the lower limits of the system's three dimensions are all taken as zero
    chargedParticleForceUnitTestProblem.setSystemSize(12.0,12.0,12.0); //(x dimension,y dimension, z dimension)

    //******************************************************ii) Particle details and geometry******************************************************
    //setting the radius of the *individual spherical particles* which will form the basis of the larger, conglomerate particles
    chargedParticleForceUnitTestProblem.setSubParticleRadius(0.25);
    //setting the total number of full *conglomerate* (hexagonal!) particles within the system
    chargedParticleForceUnitTestProblem.setParticleNumber(3);

    //setting the number of particles in the hexagon's initial row
    //can be adjusted to create particles of subtly differing shapes and elongations
    chargedParticleForceUnitTestProblem.setRowSize(2);
    //the total number of rows to include in the hexagon
    //Note: to correctly form a hexagonal particle, should always be an odd integer
    chargedParticleForceUnitTestProblem.setRowNumber(7);


    //setting the position within the system of the first conglomerate particle
    //note that the value entered corresponds to the placement top-leftmost sphere in the hexagon
    chargedParticleForceUnitTestProblem.setBasePosition(Vec3D (2.0,0.5,6.0));
    //setting the displacement by which all subsequent particles are placed
    chargedParticleForceUnitTestProblem.setMacroSeparation(Vec3D(2.0, 3.0, 0.0));

    //*********************************************************************************************************************************************
    //********************************************************Setting Interaction Details**********************************************************
    //*********************************************************************************************************************************************

    //********************************************************i)For charged interactions***********************************************************
    //setting the maximal strength and range of long-range forces between charged particles
    //Note that the force range is measured from the *edge* of particles
    chargedParticleForceUnitTestProblem.setChargeForceAndRange(3.0,4.0); //(force,range)
    //Putting details in more user-friendly terms
    //the maximum force exerted, i.e. when particles are in contact
    double maximumForce = 3.0;
    //the range of the force.
    //Note that the range is measured from the EDGE of the particle!
    double forceRange = 4.0;
    //based on user inputs, calculated the necessary adhesion stiffness
    double adStiffness = maximumForce / forceRange;
    //setting the charge values for individual particle species
    //Note that value must either be 1 (+ve charge), -1 (-ve charge)
    //or zero (no charge)
    //the charge for the first species...
    double chargeOne = 1;
    //...and that of the second species
    double chargeTwo = -1;
    //setting parameters for the closer-range van der Waals forces
    //Note that, for realistic effect, 'waalsForceMax' must be
    //*larger than the maximal charged force* (for repulsive interactions)
    double waalsForceMax = 0;
    double waalsForceRange = 0.5;
    double waalsStiffness = waalsForceMax / waalsForceRange;

    //********************************************************ii)For particle-gluing interactions*****************************************************
    //the strength of the force holding particles together
    //(Affects separation of 'glued' particles)
    double glueStrength = 250.0;
    //The dissipation associated with glued particle interactions. Can be used to damp vibration within conglomerate particles.
    double glueDiss = 0.2;

    //********************************************************iii)For collisional interactions********************************************************
    //The stiffness of the particles themselves
    double stiffness = 10000 * adStiffness;
    //The friction coefficient for inter-particle interactions
    //Note that, for simplicity, the sliding and rolling coefficients for all particles and all walls are set equal as inter-particle contact
    //interactions are not important for the systems explored here.
    //This can, however, easily be changed if desired.
    double mu = 1.0;

    //************************************************************************************************************************************************
    //********************************************************ASSIGNING GENERAL SIMULATION PARAMETERS*************************************************
    //************************************************************************************************************************************************
    //Giving a name for the output file
    chargedParticleForceUnitTestProblem.setName("ClayParticles");
    //this parameter can be used to manually adjust the time step in order to find a balance between anccuracy and efficiency
    //NOTE: should never be higher than unity!
    //This may be necessary to correctly resolve the high stiffness and short range associated with the van der Waals force.
    double timeStepReduction = 1; //e.g. 0.1 will reduce the time step tenfold!
    //Setting the desired gravitational acceleration (or lack thereof!)
    chargedParticleForceUnitTestProblem.setGravity(Vec3D(0,0,-0));
    //setting the duration of the simulation in "simulation seconds" (determined by
    //the dimensions used in setup)
    chargedParticleForceUnitTestProblem.setTimeMax(30.);

    //************************************************************************************************************************************************
    //********************************************************SETTING DONE - YOU CAN IGNORE THE REST!!************************************************
    //************************************************************************************************************************************************
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
    species1->setSlidingFrictionCoefficient(mu);
    species1->setSlidingStiffness(2.0/7.0*stiffness);
    species1->setRollingFrictionCoefficient(mu);
    species1->setRollingStiffness(2.0/5.0*stiffness);
    //setting the interaction force properties
    //setting adhesion stiffness and force max as required to give the
    //user-requested force range.
    species1->setAdhesionForceMax(maximumForce);
    species1->setAdhesionStiffness(adStiffness);
    //Setting the charge of a particle:
    //Positive (1) negative (-1) or neutral (0)
    species1->setCharge(chargeOne);
    //Assigning bond properties:
    //The maximum force exerted by the 'bond'
    //(determines the size of the overlap)
    species1->setBondForceMax(glueStrength);
    //the dissipation of the interaction -
    //used to damp out any unphysical oscillations in the contacting pair
    species1->setBondDissipation(glueDiss);
    //setting van der Waals force
    species1->setVanDerWaalsForceMax(waalsForceMax);
    species1->setVanDerWaalsStiffness(waalsStiffness);


    //SPECIES 2
    //setting the material properties of the particles
    species2->setDensity(6./constants::pi);
    species2->setStiffness(stiffness);
    species2->setSlidingFrictionCoefficient(mu);
    species2->setSlidingStiffness(2.0/7.0*stiffness);
    species2->setRollingFrictionCoefficient(mu);
    species2->setRollingStiffness(2.0/5.0*stiffness);
    //setting the interaction force properties
    //setting adhesion stiffness and force max as required to give the
    //user-requested force range
    species2->setAdhesionForceMax(maximumForce);
    species2->setAdhesionStiffness(adStiffness);
    //Setting the charge of a particle:
    //Positive (1) negative (-1) or neutral (0)
    species2->setCharge(chargeTwo);
    //Assigning bond properties:
    //The maximum force exerted by the 'bond'
    //(determines the size of the overlap)
    species2->setBondForceMax(glueStrength);
    //the dissipation of the interaction -
    //used to damp out any unphysical oscillations in the contacting pair
    species2->setBondDissipation(glueDiss);
    species2->setVanDerWaalsForceMax(waalsForceMax);
    species2->setVanDerWaalsStiffness(waalsStiffness);

    //FOR MIXED (SPECIES1-SPECIES2) INTERACTIONS
    //adding also the ability to alter the mixed-particle-interaction properties
    auto species1_2 = dynamic_cast<LinearViscoelasticFrictionChargedBondedMixedSpecies*>(chargedParticleForceUnitTestProblem.speciesHandler.getMixedObject(species1->getIndex(),species2->getIndex()));
    //for now, simply setting mixed object parameters equal to those of species 1
    //setting the material properties to the average values of the 2 particle species undergoing interaction
    species1_2->mixAll(species1,species2);
    species1_2->setSlidingFrictionCoefficient(mu);
    species1_2->setSlidingStiffness(2.0/7.0*stiffness);
    species1_2->setRollingFrictionCoefficient(mu);
    species1_2->setRollingStiffness(2.0/5.0*stiffness);
    species1_2->setVanDerWaalsForceMax(waalsForceMax);
    species1_2->setVanDerWaalsStiffness(waalsStiffness);

    //*************************************ASSIGNING PROPERTIES TO WALL SPECIES*********************************************
    //SPECIES W
    //setting the material properties of the particles
    speciesW->setDensity(6./constants::pi);
    speciesW->setStiffness(stiffness);
    speciesW->setSlidingFrictionCoefficient(mu);
    speciesW->setSlidingStiffness(2.0/7.0*stiffness);
    speciesW->setRollingFrictionCoefficient(mu);
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
    speciesW->setVanDerWaalsForceMax(0.0);
    speciesW->setVanDerWaalsStiffness(0.0);


    //FOR MIXED (PARTICLE-WALL) INTERACTIONS
    //walls and species 1 particles
    auto speciesW_1 = dynamic_cast<LinearViscoelasticFrictionChargedBondedMixedSpecies*>(chargedParticleForceUnitTestProblem.speciesHandler.getMixedObject(speciesW->getIndex(),species1->getIndex()));
    //setting the material properties to the average values of the 2 particle species undergoing interaction
    speciesW_1->mixAll(speciesW,species1);
    speciesW_1->setSlidingFrictionCoefficient(mu);
    speciesW_1->setSlidingStiffness(2.0/7.0*stiffness);
    speciesW_1->setRollingFrictionCoefficient(mu);
    speciesW_1->setRollingStiffness(2.0/5.0*stiffness);
    speciesW_1->setVanDerWaalsForceMax(0.0);
    speciesW_1->setVanDerWaalsStiffness(0.0);


    //walls and species 1 particles
    auto speciesW_2 = dynamic_cast<LinearViscoelasticFrictionChargedBondedMixedSpecies*>(chargedParticleForceUnitTestProblem.speciesHandler.getMixedObject(speciesW->getIndex(),species2->getIndex()));
    //setting the material properties to the average values of the 2 particle species undergoing interaction
    speciesW_2->mixAll(speciesW,species2);
    speciesW_2->setSlidingFrictionCoefficient(mu);
    speciesW_2->setSlidingStiffness(2.0/7.0*stiffness);
    speciesW_2->setRollingFrictionCoefficient(mu);
    speciesW_2->setRollingStiffness(2.0/5.0*stiffness);
    speciesW_2->setVanDerWaalsForceMax(0.0);
    speciesW_2->setVanDerWaalsStiffness(0.0);

    //Calculating particle mass and collision time (tc) from the parameters manually provided by the user.
    double mass = species1->getMassFromRadius(chargedParticleForceUnitTestProblem.getSubParticleRadius());
    double tc = species1->getCollisionTime(mass);
    chargedParticleForceUnitTestProblem.setTimeStep(0.02*tc*timeStepReduction);

    //setting parameters for visualisation using the xballs software package
    chargedParticleForceUnitTestProblem.setXBallsAdditionalArguments("-cmode 8 -solidf -v0");
    //solving the problem!
    chargedParticleForceUnitTestProblem.solve(argc,argv);

    std::vector<BaseParticle*>::iterator pIt = chargedParticleForceUnitTestProblem.particleHandler.begin();

}

























