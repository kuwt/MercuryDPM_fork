//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Logger.h"

/// In this file, the rolling behaviour of the tangential spring is tested. This is done by placing one normal partilce on top of a fixed partilce and letting graviry roll it over the other particle until it loses contact.
class RodsEF_2D : public DPMBase {
public:
	
	Mdouble newX, newY, newZ;
	Mdouble epsilonX, epsilonY, epsilonZ; 
	
    void setupInitialConditions()
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
        BaseParticle P0,P1,P2,P3,P4,P5;

        //setting the initial overlap of bonded particles forming a single
        //conglomerate particle
        auto species = dynamic_cast<LinearViscoelasticFrictionChargedBondedSpecies*>(speciesHandler.getObject(0));
        double delta = species->getBondForceMax()/species->getStiffness();
        //Vec3D microSeparationInitial(radius+radius-delta,0.0,0.0);
        double microSeparationInitialx=radius+radius-delta;
        double microSeparationInitialy=0.0;
        double microSeparationInitialz=0.0;

        //Vec3D microSeparation(0.0,0.0,0.0);
        double microSeparationx=0.0;
        double microSeparationy=0.0;
        double microSeparationz=0.0;

        //setting the overlap between rows of particles forming a hexagon
        double microOverlap = 0.02;

        //creating a counter variable to keep track of the objects being created
        int ob = 0;

        //****************************************For simple rods**********************************************************************************
	//Setting up a single, falling rod
        
            P0.setSpecies(speciesHandler.getObject(0));
            P1.setSpecies(speciesHandler.getObject(1));
            P2.setSpecies(speciesHandler.getObject(1));
            P3.setSpecies(speciesHandler.getObject(1));
            P4.setSpecies(speciesHandler.getObject(0));
	    

            Vec3D microSeparation;

	    //Setting first particle position
	    double xx= (getXMin() + getXMax()) / 2.0;
	    double yy= (getYMin() + getYMax()) / 2.0;
	    double zz= (getZMax() - 2.0);

            double Vx= 0.0;
	    double Vy= 0.0;
	    double Vz= 0.0;
			
            P0.setPosition(Vec3D(xx + microSeparationx,yy,zz + microSeparationz));
            microSeparationx += microSeparationInitialx;
            microSeparationz += microSeparationInitialz;
            
            P1.setPosition(Vec3D(xx + microSeparationx,yy,zz + microSeparationz));
            microSeparationx += microSeparationInitialx;
            microSeparationz += microSeparationInitialz;
            
            P2.setPosition(Vec3D(xx + microSeparationx,yy,zz + microSeparationz));
            microSeparationx += microSeparationInitialx;
            microSeparationz += microSeparationInitialz;
            
            P3.setPosition(Vec3D(xx + microSeparationx,yy,zz + microSeparationz));
            microSeparationx += microSeparationInitialx;
            microSeparationz += microSeparationInitialz;
            
            P4.setPosition(Vec3D(xx + microSeparationx,yy,zz + microSeparationz));

            //setting initial velocity, if desired
	    /* WHY WAS THIS ORIGINALLY ANGULAR??
            P0.setAngularVelocity(Vec3D(Vx,Vy,Vz));
            P1.setAngularVelocity(Vec3D(Vx,Vy,Vz));
            P2.setAngularVelocity(Vec3D(Vx,Vy,Vz));
            P3.setAngularVelocity(Vec3D(Vx,Vy,Vz));
            P4.setAngularVelocity(Vec3D(Vx,Vy,Vz));
	    */
	    P0.setVelocity(Vec3D(Vx,Vy,Vz));
            P1.setVelocity(Vec3D(Vx,Vy,Vz));
            P2.setVelocity(Vec3D(Vx,Vy,Vz));
            P3.setVelocity(Vec3D(Vx,Vy,Vz));
            P4.setVelocity(Vec3D(Vx,Vy,Vz));


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
        

	    //Creating a "blade" upon which falling rod will hopefully shatter
	    P5.setSpecies(speciesHandler.getObject(0));
	    //P5.setPosition(Vec3D(((getXMax() + getXMin())/2.0)+1,yy,2.0));
	    P5.setPosition(Vec3D(((getXMax() + getXMin())/2.0)+1.5,yy,2.0));
	    P5.setVelocity(Vec3D(0,0,0));
	    P5.setRadius(0.125);
	    P5.fixParticle();
	    particleHandler.copyAndAddObject(P5);
        //***************************************Creating Walls and Assigning Species and Properties**********************************************

        //*******************************************************For Periodic Walls***************************************************************

        if (!solidWalls) {
            PeriodicBoundary zBounds;
            zBounds.set(Vec3D(0.0, 0.0, 1.0), getZMin(), getZMax());
            boundaryHandler.copyAndAddObject(zBounds);

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
    int nParticles = 1;

    //the position of the first particle put into the system
    Vec3D basePosition;

    //a parameter determining how particles are spread throughout the system
    Vec3D macroSeparation;

    //the size of the 'first row' of a hexagon being created
    int nWidthInitial;
    //The total number of rows forming the hexagon
    int nRow;


};

Logger<Log::ERROR> testLogger("RodsEF_2D");

int main(int argc, char *argv[])
{

    //******************************************************i) System walls and geometry***********************************************************
    
    RodsEF_2D problem;
  
    problem.setSolidWalls(true);
    problem.setSystemSize(10.0,10.0,10.0); 

    //******************************************************ii) Particle details and geometry******************************************************
    
    problem.setSubParticleRadius(0.5);
   
    //********************************************************Setting Interaction Details**********************************************************
    
    //********************************************************i)For charged interactions***********************************************************
    
    problem.setChargeForceAndRange(10.0,1.0); 
    double maximumForce = 10.0;
    double forceRange = 1.0;
    double adStiffness = maximumForce / forceRange;
    double chargeOne = 0;
    double chargeTwo = 0;
    double waalsForceMax = 1e-9;
    double waalsForceRange = 1e-9;
    double waalsStiffness = waalsForceMax / waalsForceRange;

    //********************************************************ii)For particle-gluing interactions*****************************************************
    //the strength of the force holding particles together
    //(Affects separation of 'glued' particles)
    double glueStrength = 0.1e4;
    //The dissipation associated with glued particle interactions. Can be used to damp vibration within conglomerate particles.
    double glueDiss = 600;

    //********************************************************iii)For collisional interactions********************************************************
    //The stiffness of the particles themselves
    double stiffness = 1.0e5;
    //NEW specifically defining the Rolling stiffness of the particles
    double rollingStiffness = 5.0e1;
    //NEW adding in torsion friction, stiffness and dissipation
    double torsionStiffness = 2.5e1;
    double mu = 1.0e-1;
    //NEW - individually specifying rolling friction
    double mu_r = 1.0e-2;
    //NEW adding in torsion friction, stiffness and dissipation
    double mu_t = 1.0e-5;
    //NEW - allowing sliding and rolling dissipation to be set
    double gamma_s = 0;
    double gamma_r = gamma_s;
    //NEW adding in torsion friction, stiffness and dissipation
    double gamma_t = gamma_r;
    //************************************************************************************************************************************************
    //********************************************************ASSIGNING GENERAL SIMULATION PARAMETERS*************************************************
    //************************************************************************************************************************************************
    problem.setName("RodsEF_2D");
    
    double timestepReduction = 0.1; //e.g. 0.1 will reduce the time step tenfold!
    problem.setGravity(Vec3D(0.0,0.0,-90));
    problem.setTimeMax(100.0);
    problem.epsilonX = 1e-1;
    problem.epsilonY = 1e-1;
    problem.epsilonZ = 1e-1;

    //************************************************************************************************************************************************
    //********************************************************SETTING DONE - YOU CAN IGNORE THE REST!!************************************************
    //************************************************************************************************************************************************
    //setting the species of particles
    auto species1 = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionChargedBondedSpecies());
    //adding second species to allow for attractive interactions
    auto species2 = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionChargedBondedSpecies());
    //setting species for system walls
    auto speciesW = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionChargedBondedSpecies());

    //*************************************ASSIGNING PROPERTIES TO PARTICLE SPECIES*****************************************
    //SPECIES 1
    //setting the material properties of the particles
    species1->setDensity(6./constants::pi);
    species1->setStiffness(stiffness);
    species1->setSlidingFrictionCoefficient(mu);
    species1->setSlidingStiffness(2.0/7.0*stiffness);
    //new
    species1->setRollingFrictionCoefficient(mu_r);
    //new
    species1->setTorsionFrictionCoefficient(mu_t);
    //species1->setRollingStiffness(2.0/5.0*stiffness);
    //New: adding individual rolling stiffness
    species1->setRollingStiffness(rollingStiffness);
    //...and torsion stiffness
    species1->setTorsionStiffness(torsionStiffness);
    //New:Adding a sliding and rolling dissipation
    species1->setSlidingDissipation(gamma_s);
    species1->setRollingDissipation(gamma_r);
    //...and torsional dissipation
    species1->setRollingDissipation(gamma_t);
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
	species1->setDissipation(50);

    //SPECIES 2
    //setting the material properties of the particles
    species2->setDensity(6./constants::pi);
    species2->setStiffness(stiffness);
    species2->setSlidingFrictionCoefficient(mu);
    species2->setSlidingStiffness(2.0/7.0*stiffness);
    //New
    species2->setRollingFrictionCoefficient(mu_r);
    //species2->setRollingStiffness(2.0/5.0*stiffness);
    //New: adding individual rolling stiffness
    //New:Adding a sliding and rolling dissipation
    species2->setSlidingDissipation(gamma_s);
    species2->setRollingDissipation(gamma_r);

    species2->setRollingStiffness(rollingStiffness);
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
    species2->setDissipation(50);

    //FOR MIXED (SPECIES1-SPECIES2) INTERACTIONS
    //adding also the ability to alter the mixed-particle-interaction properties
    auto species1_2 = dynamic_cast<LinearViscoelasticFrictionChargedBondedMixedSpecies*>(problem.speciesHandler.getMixedObject(species1->getIndex(),species2->getIndex()));
    //for now, simply setting mixed object parameters equal to those of species 1
    //setting the material properties to the average values of the 2 particle species undergoing interaction
    species1_2->mixAll(species1,species2);
    species1_2->setSlidingFrictionCoefficient(mu);
    //species1_2->setSlidingStiffness(2.0/7.0*stiffness);
    species1_2->setSlidingStiffness(stiffness);
    //NEW
    species1_2->setRollingFrictionCoefficient(mu_r);
    //species1_2->setRollingStiffness(2.0/5.0*stiffness);
    //NEW
    species1_2->setRollingStiffness(rollingStiffness);
    //New:Adding a sliding and rolling dissipation
    species1_2->setSlidingDissipation(gamma_s);
    species1_2->setRollingDissipation(gamma_r);

    species1_2->setVanDerWaalsForceMax(waalsForceMax);
    species1_2->setVanDerWaalsStiffness(waalsStiffness);
    //species1_2->setDissipation(0.5);

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
    speciesW->setDissipation(50);


    //FOR MIXED (PARTICLE-WALL) INTERACTIONS
    //walls and species 1 particles
    auto speciesW_1 = dynamic_cast<LinearViscoelasticFrictionChargedBondedMixedSpecies*>(problem.speciesHandler.getMixedObject(speciesW->getIndex(),species1->getIndex()));
    //setting the material properties to the average values of the 2 particle species undergoing interaction
    speciesW_1->mixAll(speciesW,species1);
    speciesW_1->setSlidingFrictionCoefficient(mu);
    speciesW_1->setSlidingStiffness(2.0/7.0*stiffness);
    speciesW_1->setRollingFrictionCoefficient(mu);
    speciesW_1->setRollingStiffness(2.0/5.0*stiffness);
    speciesW_1->setVanDerWaalsForceMax(0.0);
    speciesW_1->setVanDerWaalsStiffness(0.0);
    //speciesW_1->setDissipation(0.5);


    //walls and species 1 particles
    auto speciesW_2 = dynamic_cast<LinearViscoelasticFrictionChargedBondedMixedSpecies*>(problem.speciesHandler.getMixedObject(speciesW->getIndex(),species2->getIndex()));
    //setting the material properties to the average values of the 2 particle species undergoing interaction
    speciesW_2->mixAll(speciesW,species2);
    speciesW_2->setSlidingFrictionCoefficient(mu);
    speciesW_2->setSlidingStiffness(2.0/7.0*stiffness);
    speciesW_2->setRollingFrictionCoefficient(mu);
    speciesW_2->setRollingStiffness(2.0/5.0*stiffness);
    speciesW_2->setVanDerWaalsForceMax(0.0);
    speciesW_2->setVanDerWaalsStiffness(0.0);
    //speciesW_2->setDissipation(0.5);

    //Calculating particle mass and collision time (tc) from the parameters manually provided by the user.
    double mass = species1->getMassFromRadius(problem.getSubParticleRadius());
    double tc = species1->getCollisionTime(mass);
    problem.setTimeStep(0.02*tc*timestepReduction);
    problem.setSaveCount(100);

    //setting parameters for visualisation using the xballs software package
    problem.setXBallsAdditionalArguments("-cmode 8 -solidf -v0");
    //solving the problem!
    problem.solve(argc,argv);

    std::vector<BaseParticle*>::iterator pIt = problem.particleHandler.begin();

}


























