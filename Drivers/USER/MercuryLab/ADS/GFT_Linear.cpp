//Copyright (c) 2013-2017, The MercuryDPM Developers Team. All rights reserved.
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
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Species/Species.h"
#include "Species/HertzianViscoelasticMindlinSpecies.h"
#include <string>
#include <sstream>

/**
 * Simulates Granular Flow Target
 */
class GFT : public Mercury3D {

  void setupInitialConditions() override {
    //setting system size
    setMin({-width, -width, 0});
    setMax({width, width, height});

    //add walls

    //add outer cylindrical wall
    AxisymmetricIntersectionOfWalls w0({0,0,0}, {0,0,1}, {}, speciesHandler.getObject(0));
    w0.addObject({1,0,0},{width,0,0});
    wallHandler.copyAndAddObject(w0);

    //add outer funnel wall
    AxisymmetricIntersectionOfWalls w1({0,0,0}, {0,0,1}, {}, speciesHandler.getObject(0));
    //Enforcing a specific funnel angle
    w1.addObject({50.0/45.0,0,-1},{fwidth,0,0});
    w1.addObject({0,0,1},{fwidth,0,0});
    wallHandler.copyAndAddObject(w1);

    //add inner wall
    AxisymmetricIntersectionOfWalls w2({0,0,0}, {0,0,1}, {}, speciesHandler.getObject(0));
    w2.addObject({-1,0,0},{0.5*width,0,iheight});
    w2.addObject({1,0,0.001},{0.5*width,0,iheight});
    wallHandler.copyAndAddObject(w2);

    // uncomment to add flat wall at base, thus preventing outflow
    // InfiniteWall w({0,0,-1}, {0,0,0}, speciesHandler.getObject(0));
    // wallHandler.copyAndAddObject(w);

    //add particles
    BaseParticle p;
    p.setSpecies(speciesHandler.getObject(0));
    //creating a random number whose magnitude is defined by 'polydispersity'
    //to give particles within the system a randomised distribution of sizes
    //within the desired bounds
    p.setRadius(radius);
    Mdouble zMax = 0;
    Mdouble dz = 0.0005*radius; //this number can be adjusted to make particle insertion quicker
    while (particleHandler.getNumberOfObjects() < particleNumber) {
      // define the position at which the particle is to be inserted
      Mdouble x, y, rr, z = random.getRandomNumber(0, zMax);
      do {
	x = random.getRandomNumber(-width, width);
	y = random.getRandomNumber(-width, width);
	rr = x * x + y * y;
      } while (rr > mathsFunc::square(width - radius) || rr < mathsFunc::square(0.5 * width + radius));
      p.setPosition({x, y, z});
      // check if particle can be inserted
      if (checkParticleForInteraction(p)) {
	particleHandler.copyAndAddObject(p);
	double pol = random.getRandomNumber(-polydispersity,polydispersity);
	p.setRadius(radius*(1+pol));
	std::cout << '+';
	std::cout.flush();
      } else {
	//gradually increase the insertion domain to to insert particles as low as possible
	zMax += dz;
	//std::cerr << "warning: failed insertion attempt" << std::endl;
      }
    }
    std::cout << " all particles inserted below z/h=" << zMax/height << std::endl;
  }

  void printTime () const override {
    std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
              << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
              << ", EneKin=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
              << std::endl;
    std::cout.flush();
  }

  ///creating a reinsertion boundary
  void actionsAfterTimeStep() override {
    //take all particles below height z=-50mm and reinsert them into the GFT above height zMax
    for (auto p: particleHandler) {
      if (p->getPosition().Z < -50e-3) {
	p->setVelocity({0, 0, 0});
	Mdouble zMax = 0;
	do {
	  zMax += 1e-3;
	  Mdouble x, y, rr, z = random.getRandomNumber(getZMax(), getZMax() + zMax);
	  do {
	    x = random.getRandomNumber(getXMin(), getXMax());
	    y = random.getRandomNumber(getYMin(), getYMax());
	    rr = x * x + y * y;
	  } while (rr > mathsFunc::square(getXMax() - p->getRadius()) ||
		   rr < mathsFunc::square(0.5 * getXMax() + p->getRadius()));
	  p->setPosition({x, y, z});
	} while (!checkParticleForInteraction(*p));
	//std::cout << '^'; //uncomment this line to get output every time a particle is reinserted
      }
    }
  }

private:

  ///a variable which 'accelerates' simulations by making the timestep larger
  ///useful for trial runs
  ///(initialized to default 1)
  Mdouble speedUp=1;
  ///The number of particles within the system
  unsigned particleNumber;
  ///the factor by which the particle size is multiplied in test simulations (set to 1 for real simulations)
  ///(initialized to default 1)
  Mdouble scaleUp=1;
  ///the radius possessed by particles
  Mdouble radius = 1.5e-3 * scaleUp;
  ///the (outer) width of the system
  Mdouble width;
  ///the funnel width
  Mdouble fwidth;
  ///the height of the inner wall
  Mdouble iheight;
  ///The height of the system as a whole
  Mdouble height;
  ///A baseline value of the number of particles which can be used to autonomously vary the number
  ///of particles in the system as the particle size scaling is altered such that a constant
  ///packing fraction is maintained.
  ///The value is created to be consistent with the values assigned in the original version of
  ///the code (N = 2800, scaleUp = 4, radius = 1.5e-3) to ensure consistency
  unsigned basalParticleNumber = 2800;
  ///The (fractional) degree of polydispersity within the system
  double polydispersity;

  //******************************get and set functions******************************************************

public:

  //setting the relevant system dimensions, namely:
  //the positioning of the *outer* cylindrical wall and
  //the *inner* cylindrical wall or 'funnel'
  void setWidths(Mdouble wo, Mdouble wi) {
    width = wo;
    fwidth = wi;
  }

  //the total height of the system and the height of the inner tube
  void setHeights(Mdouble ht, Mdouble hi) {
    height = ht;
    iheight = hi;
  }

  //functions to 'get' useful system size parameters
  double getWidth() {
    return width;
  }
  double getInnerWidth() {
    return fwidth;
  }
  double getHeight() {
    return height;
  }
  double getInnerHeight() {
    return iheight;
  }
  //setting the desired number or particles
  void setParticleNumber(unsigned n) {
    particleNumber = n;
  }

  //...and making it possible to retrieve particle numbers
  double getParticleNumber() {
    return particleNumber;
  }

  //a function to set the desired 'speedUp', i.e. the increase
  //in timestep used to make trial simulations faster...
  void setSpeedUp(double s) {
    speedUp = s;
  }

  //...and a function to return the value of the speed-up
  double getSpeedUp() {
    return speedUp;
  }

  //a function to set the particle radius and scaling factor
  //(setting both will make it easier to create a scaling for the
  //total system particle number)
  void setRadiusAndScaling(double r, double s) {
    radius = r * s;
    scaleUp = s;
    //adjusting also the number of particles to ensure an equal packing
    //fraction as the scaling is adjusted
    //the adjustment is based on the original parameters with which the
    //simulation was created - r = 1.5e-3, N = 2800, scaleUp = 4
    particleNumber = basalParticleNumber * (0.0015) * 4.0 * (0.0015) * 4.0 * (0.0015) * 4.0
      / (radius * radius * radius  );
    //note that 'basalParticleNumber' is left as an external variable so that
    //the packing of the system can still be varied if desired

    std::cout << "Requested number of particles:" << particleNumber << std::endl;
  }

  //a get function to allow radius to be used in calculations elsewhere
  double getRadius() {
    return radius;
  }

  //a get function for the scale-up used on particles
  double getScaleUp() {
    return scaleUp;
  }

  //a function to set the degree of polydispersity within the system
  void setPolydispersity( double poly ) {
    polydispersity = poly;
  }
  //a function to retrieve the value of the polydispersity of particles
  double getPolydispersity() {
    return polydispersity;
  }
};

int main(int argc UNUSED, char *argv[] UNUSED) {
  //setting up the problem
  GFT gft;
  std::string nameRoot = "hertzianAngleTest"; //setting up the 'root' of the files name
  //*********************************************************************************************************
  //***********************************Major User-Defined Variables******************************************
  //*********************************************************************************************************

  //*********************************Physical properties of particles****************************************
  HertzianViscoelasticMindlinSpecies species;
  Mdouble elasticModulus = 4.11e11;
  Mdouble poissonRatio = 0.28;
  species.setDensity(19250); //implementing accurate value of particle density
  species.setElasticModulusAndRestitutionCoefficient(elasticModulus,
						     0.95); //(Elastic modulus, restitution coefficient) (normally 0.1)
  species.setSlidingFrictionCoefficient(0.225);//setting the sliding coefficient of friction (normally 0.5)
  species.setShearModulus( (0.5 * 2.0 * elasticModulus / (1 + poissonRatio) ) / 2.0); //setting the shear modulus (weird numbers highlight that we are using the effective value)
  species.setSlidingDissipation(0.5);
  gft.speciesHandler.copyAndAddObject(species);
  gft.setRadiusAndScaling(1.5e-3, 1.66666666); //(particle radius,scale factor)
  //gft.setParticleNumber(5); // uncomment if you want to set the particle number manually
  gft.setPolydispersity(0.0); // setting the fractional polydispersity of the system's particles
  //*************************************System Size and Details*********************************************
  gft.setWidths(60e-3, 15e-3); //(outer cylinder radius, outlet radius) Note inner cylinder is 1/2 outer cylinder width
  gft.setHeights(500e-3, 100e-3); //(total system height,funnel height) **USED TO SET HOPPER ANGLE**
  gft.setGravity(Vec3D(0, 0, -9.8));//setting gravity


  //***************************************Simulation Details***********************************************
  gft.setTimeMax(3.0); //the total duration of the simulation
  gft.setSaveCount(200); //setting how often to output data to file
  //gft.fStatFile.setSaveCount(10); //setting how often to output fstat data to file
  Mdouble maxRelativeVelocity = 1e-2; //used to accurately calculate contact time
  gft.setSpeedUp(3.0); //setting the required speed-up for use when testing code (set 1 for real simulations!)
  Mdouble tc = species.getCollisionTime(2.0 * gft.getRadius(), species.getDensity(),
					maxRelativeVelocity ); //setting contact time based on particle properties
  logger(INFO, "Collision time %", tc);
  gft.setTimeStep(gft.getSpeedUp() * 0.02 * tc);//setting the timestep of the simulation, dependent on the desired 'speedUp'

  //*********************************************************************************************************
  //*************************************************END*****************************************************
  //*********************************************************************************************************

  //****************************Dependent (not directly user-defined) parameters*****************************
  //creating an informative name for the output file
  //creating an empty stringstream to read in relevant information
  std::stringstream fullName(nameRoot);
  //adding all relevant information to the stringstream
  fullName << nameRoot << "-scaleUp" << gft.getScaleUp() << "-speedUp" << gft.getSpeedUp()
	   << "-N" << gft.getParticleNumber() << "-mu" << species.getSlidingFrictionCoefficient()
	   << "-poly" << gft.getPolydispersity()
	   << "-outerW" << gft.getWidth() << "-innerW"  << gft.getInnerWidth()
	   << "-heightTot" << gft.getHeight() << "-heightIn" << gft.getInnerHeight() << "new";
  //naming the output file
  gft.setName(fullName.str());
  gft.setXBallsAdditionalArguments("-v0 -solidf");
  //running the simulation!
  gft.solve();
  helpers::writeToFile("GFT.gnu",
		       "p [][] 'GFT.fstat' u 8:10 w p");
  logger(INFO, "type 'gnuplot GFT.gnu' to check force law");
  return 0;
}
