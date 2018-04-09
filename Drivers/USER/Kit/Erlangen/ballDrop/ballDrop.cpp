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

//based on /storage2/usr/people/sluding/MDCC/C3DshearXL30/MU0_LONG2
#include "Mercury3D.h"
#include "Species/HertzianViscoelasticSlidingFrictionSpecies.h"
#include "Particles/BaseParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"


/** This self test was written to test the speed of particle creation in MercuryDPM.
 * The code below creates 200'000 non-overlapping particles of homogeneous size.
 * ///\todo is the above still true? it seems that there are 20^3 particles are constructed
 * ///\todo this code still feels a bit messy and unclear
 **/
class BallDrop : public Mercury3D
{

private:
   
    //The radius of the large particle
    Mdouble bigBallRadius;

    //Boolean to determine if the large particle has been added
    bool bigBallAdded = false;

    //The initial position within the system at which an attempt is made to insert the large particle
    Mdouble tryHeight = 0;
    

    //Threshold system KE below which the system can be assumed 'settled'
    // Makes sense to take this as a value per-particle to remove any system
    //size dependence
    double thresholdKE = 1e-5;
 
    //The time at which to begin checking if KE is decayed
    double checkTime = 1.0;

    //The interval at which intial KE checks are performed
    double checkTimeIncrement = 0.5;

    // The dimensions of the 'bed', i.e. the backfill area
    Mdouble bedWidth = 20.0; //length of backfill area (x)
    Mdouble bedDepth = 20.0; //depth of backfill area (y)
    Mdouble bedHeight = 10.0; // height of confining walls (z)

    //allowing the user to determine the relative size & position of the area from which
    //particles are inserted into the system
    //Note these values are all relative to the bed itself.
    Mdouble fill_xMin = 0;
    Mdouble fill_xMax = bedWidth;
    Mdouble fill_yMin = 0;
    Mdouble fill_yMax = bedDepth;
    Mdouble fill_zMin = bedHeight;
    Mdouble fill_zMax = 2*bedHeight;

    Mdouble vRatio = 1.0; //the ratio of inlet volume to system volume 
			  //default 1, i.e. inlet volume = system volume

    //A boolean to allow the user to choose whether they want a "full" bed or whether
    //they want to choose the number of particles in the system
    bool fullSystem = true;

    //A double which determines the relative size of the fill volume to the actual intended system volume
    double overFill = 1.0;

    //The number of particles to include if the user chooses not to fill the system
    unsigned int N = 100;

    //The radius of a particle 
    Mdouble rad = 0.5;
    //the radius of the smallest possible particle
    Mdouble radMin;
    //the particle density
    Mdouble rho = 6.0 / constants::pi;
    //the particle stiffness
    Mdouble k = 500;
    //the restitution coefficient
    Mdouble eps = 0.1;


    //The gravitational acceleration
    Mdouble g = -1.0;

    //Frictional properties
    Mdouble mu_s = 0.5; //sliding friction
    Mdouble mu_r = mu_s / 10.0; //rolling friction

    //TO DO - put in sanity check for particles 'missing' the bed!
public:
    
    //constructor
    BallDrop(Mdouble W, Mdouble D,Mdouble H, Mdouble r, Mdouble sizeDistribution, int particleNumber, Mdouble grav,
	     Mdouble slidingFriction, Mdouble rollingFriction, Mdouble timeStop, 
	     Mdouble density, Mdouble stiffness, Mdouble restitution, Mdouble OF, Mdouble bigRad) :
            sizeDistribution_(sizeDistribution)
    {

	bigBallRadius = bigRad;

	bedWidth = W;
	bedDepth = D;
	bedHeight = H;

        setXMin(0);
        setXMax(bedWidth);
        setYMin(0);
        setYMax(bedDepth);
        setZMin(0);
        setZMax(bedHeight);

	//Setting the initial height at which particle insertion is attempted to the bottom of the system 
	// plus one particle radius
	tryHeight = getZMin() + bigBallRadius;

	g = grav;
        setGravity({0.0,0.0,g});

	rad = r;

	N = particleNumber;

	mu_s = slidingFriction;
	mu_r = rollingFriction;

	rho = density;
	k = stiffness;

	eps = restitution;

	overFill = OF;

        logger(INFO, "Creating system of width [%,%], depth [%,%] height [%,%]",
                        getXMin(), getXMax(), getYMin(),getYMax(), getZMin(), getZMax());

        //create new species
        species = speciesHandler.copyAndAddObject(HertzianViscoelasticSlidingFrictionSpecies());
        species->setDensity(rho);

	//setting contact properties
	//determining the smallest possible radius
	radMin = cbrt(rad / (sizeDistribution_ * sizeDistribution_ + 2*rad) / (sizeDistribution_ + 2*rad));
	
	//for safety, setting the species properties according to the smallest possible particle 
        species->setElasticModulusAndRestitutionCoefficient(1e7, eps);


	//setting frictional properties
	species->setSlidingStiffness(species->getStiffness()*2./7.);
        species->setSlidingFrictionCoefficient(mu_s);
        species->setSlidingDissipation(species->getDissipation()*2./7.);
	species->setRollingStiffness(species->getStiffness()*2.0/7.0);
        species->setRollingFrictionCoefficient(mu_r);
        species->setRollingDissipation(species->getDissipation()*2./7.);

	//for safety, setting collision time, and hence timestep, according to smallest particles present
        setTimeStep(0.02 * species->getCollisionTime(species->getMassFromRadius(radMin)));
        logger(INFO, "time step used %", getTimeStep());

        setTimeMax(timeStop); //setting the duration of the simulation

        //set walls
        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0)); // setting species for all walls
	//Adding a flat basewall
        w.set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin())); 
        wallHandler.copyAndAddObject(w);
	// Adding a curved outer wall
	Vec3D drumCenter = {0.5*(getXMin() + getXMax()), 0.5*(getYMin() + getYMax()), 0.5*(getZMin() + getZMax())};
	auto drumWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        drumWall -> setSpecies(species);
        drumWall -> setPosition(drumCenter);
        drumWall -> setOrientation(Vec3D(0.0,0.0,1.0));
        drumWall -> addObject(Vec3D(1,0,0), Vec3D(getXMax()/2.0,0.0,0.0));
    }

    void setupInitialConditions()
    {
        //number of particles to be inserted
        unsigned int n;
        //determining if to create full system or user-defined n
	if (fullSystem) {
            n = static_cast<unsigned int>((getXMax() - getXMin()) * (getYMax() - getYMin()) *
                                          (getZMax() - getZMin()) / (M_PI * 4 *  (rad * rad * rad)/ 3));
        } else {
	    n = N;
	}

        //try to find new insertable particles
        unsigned int numberOfParticlesInserted = 0; //counter recording current particle number
        BaseParticle p; //creating a particle
        p.setSpecies(species); //assigning its species
        const Mdouble s = sizeDistribution_;
        const Mdouble rMin = radMin;
        logger(INFO, "Inserting a maximum of %  particles with  %"
        "<d<%  (sizeDistribution %)", n, 2.0 * rMin, 2.0 * s * rMin, s);
        //this already changes the largest particle (before it gets inserted into DPMBase)
        p.setRadius(s * rMin);
        Vec3D position;
        unsigned int failed = 0;
        while (numberOfParticlesInserted < n)
        {
            position.X = random.getRandomNumber(getXMin() + p.getRadius(),
                                                getXMax() - p.getRadius());
            position.Y = random.getRandomNumber(getYMin() + p.getRadius(),
                                                getYMax() - p.getRadius());
            position.Z = random.getRandomNumber(getZMin() + p.getRadius(),
                                                getZMax()*overFill - p.getRadius());
            p.setPosition(position);
            if (checkParticleForInteraction(p)) //true if there is no interaction
            {
                logger(DEBUG, "insert r= %, N=%, R=%, %, %", p.getRadius() ,particleHandler.getNumberOfObjects(),
                    particleHandler.getLargestParticle()->getRadius(), particleHandler.getLargestParticle(), &p);
                particleHandler.copyAndAddObject(p);
                if (hGridNeedsRebuilding())
                {
                    hGridRebuild();
                }
                logger(DEBUG, "done inserting, N= %", particleHandler.getNumberOfObjects());
		//now resetting the particle size to a random value within distribution
                p.setRadius(random.getRandomNumber(rMin, s * rMin));
		p.setVelocity({0,0,0});
		logger(DEBUG, "next particle r= %", p.getRadius());
                numberOfParticlesInserted++;
                if (particleHandler.getNumberOfObjects() % 1000 == 0)
                {
                    logger(INFO, " %k particles inserted", particleHandler.getNumberOfObjects() / 1000);
                }
                failed = 0;
            }
            else
            {
                failed++; //if insertion fails, increments the fail counter
                if (failed > 10000) //if fails repeatedly, gives up inserting particles (assume system is full)
                {
                    logger(DEBUG, "\nfailed to insert 10000 particles in a row; aborting");
                    break;
                }
            }
        }
        logger(INFO, "\nInserted % particles", particleHandler.getNumberOfObjects());
    }


    void actionsBeforeTimeStep() {
	//Determining whether to add large particle
	//First, checking it has not already been added
	if (!bigBallAdded) {
	    //Checking if we are beyond the required minimum time for the next insertion attempt
	    if (getTime() > checkTime) {
		std::cout << "Bed relaxing. Current total kinetic energy = " << getKineticEnergy() << std::endl;
		//Checking if the system's current KE is below the threshold required to insert the particle
		if (getKineticEnergy() < thresholdKE) {
			//If so, inserting particle and alerting user
			std::cout << "Bed relaxed; inserting large particle." << std::endl;
			//Creating an infinite loop over which insertion is repeatedly attempted until
			//successful
			while(true) {
				//Optional user warning
				std::cout << "Attempting to insert particle at z = " << tryHeight << std::endl;				

				BaseParticle P; //declaring a generalised particle with no properties
				P.setSpecies(species); //assigning particle a species
				P.setRadius(bigBallRadius); //assigning chosen large particle size
				Vec3D position; //creating a vector to store particle position
				position.X = (getXMin() + getXMax()) / 2.0; //setting X position to system centre
				position.Y = (getYMin() + getYMax()) / 2.0; //setting Y position to system centre
				position.Z = tryHeight; //setting Z position to current trial value
				//Assigning current position to potential particle
				P.setPosition(position);
				//Assigning zero initial velocity
				P.setVelocity(Vec3D(0.0,0.0,0.0));
				//seeing if the particle can be inserted at the current tryHeight
				if(checkParticleForInteraction(P)) {
					//If so, inserting the particle...
					particleHandler.copyAndAddObject(P);
					//Once particle is successfully inserted, setting 'bigBallAdded' flag to true
                        		bigBallAdded = true;
					//...Letting user know the insertion height...
					std::cout << "Large particle inserted at z = " << tryHeight << std::endl;
					//...and breaking the loop
					break;
				}
				else {
					//If insertion fails at current height, incrementing the trial height by
					//one (normal-size) particle radius for next attempt...
					tryHeight = tryHeight + rad;	
				}
			}
		}
		else {
			//if not, resetting the check time to the next desired increment
			checkTime = getTime() + checkTimeIncrement; 
		}
	    }
	}
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////Non-standard functions///////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    //Function to set up the size and position of the inlet through which particles are 'poured'
    void setInlet(double w,double d, double xStart, double yStart) {
	fill_xMin = xStart; //setting start of inlet to chosen x position
	fill_xMax = xStart + w; //setting width of inlet to w
	fill_yMin = yStart; //setting start of inlet to chosen x position
        fill_yMax = yStart + d; //setting width of inlet to w
	
	//sanity check - making sure particles are not poured outside the system!
	if (fill_xMin < 0 || fill_yMin < 0 || fill_xMax > bedWidth || fill_yMax > bedDepth ) {
	    std::cerr << "Error: particle inlet lies outside of system boundaries" << std::endl;
	    exit(EXIT_FAILURE);
	}

	//ensuring the z-dimension is such that the inlet size equals vRatio times the main system size 
	//(i.e. will accurately fill the system)
	double v = (bedWidth) * (bedDepth) * (bedHeight); //volume of backfill system
	double h = v / (w * d);
	fill_zMin = bedHeight; //starting fill from (just) above bed height
	fill_zMax = bedHeight + vRatio * h;

	std::cout << "Fill region height = " << h << std::endl;
	std::cout << "Relative volume (fill region / system) = "  << w * d * h / (bedWidth * bedDepth * bedHeight) << std::endl; 
	
    }
private:
    HertzianViscoelasticSlidingFrictionSpecies* species;
    Mdouble sizeDistribution_;
};

int main(int argc, char* argv[])
{

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////System properties///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // The dimensions of the system
    const Mdouble width = 20.0; //in particle diameters
    const Mdouble depth = width; //in particle diameters
    const Mdouble height = 20.0; //in particle diameters

    const Mdouble g = -1.0; //the strength of gravity (by definition in -z direction, but can be changed)

    const Mdouble overFill = 2.0; // The factor by which the initial particle-loading volume is greater
					//than the final volume, to ensure a correct fill height

    const Mdouble N = 10; //the number of particles in the system - only used if non-full system is requested

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////Particle properties/////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    const Mdouble bigBallRadius = 2.0;
    const Mdouble r = 0.5; //the particle radius
    const Mdouble sizeDistribution = 2.0; //The factor by which particle sizes can vary
    const Mdouble mu_s = 0.5; //sliding friction
    const Mdouble mu_r = mu_s / 10.0; //rolling friction - used in linear version only
    const Mdouble rho = 6.0 / constants::pi; //particle density 
    const Mdouble k = 500; //particle stiffness
    const Mdouble e = 0.2; //restitution coefficient
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////Simulation properties///////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    const Mdouble tMax = 200.0;

    BallDrop BD(width,depth,height,r,sizeDistribution,N,g,
		mu_s,mu_r,tMax,rho,k,e,overFill,bigBallRadius);
    BD.setName("Test1");
    BD.setSaveCount(50); //records every Nth calculated timestep
    BD.setXBallsAdditionalArguments("-cmode 0 -solidf -v0");
    BD.solve(argc, argv);
    return 0;
}
