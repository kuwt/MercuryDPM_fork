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

#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Boundaries/PeriodicBoundary.h"
#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include <math.h>
#include <stdlib.h>


class Box : public Mercury3D
{
public:
	
	double r0;
	double r1;
	double B;
	double argv4;
	double argv5;
	double argv6;
	double argv7;
	double argv8;

	void setupInitialConditions()
	{
		auto species = new LinearViscoelasticSpecies; // Sets the base species for the particles and walls
	
		speciesHandler.addObject(species); // Adds the species to the setup 
		setParticleDimensions(3);
	
		double YsMod = 3.0; // Young Modulus of Nylon
		//Supposed to be 3 GPa = 3e9 Pa??

		double stiffness = (M_PI/2.0) * r0 * YsMod * 100 * 50; // Stiffness of Nylon
		// Need proper equation for this

		double density = 1150.0; // Denisty of Nylon
		double mass = M_PI * density * r0*r0*r0 * 4.0/3.0; // Mass of each bead

		species -> setDensity(density); // Sets the density of the species 
		species -> setStiffnessAndRestitutionCoefficient(stiffness, 0.9, mass); // Set stiffness and restituion coefficient of the species
		// Need true value for restitution coefficient!

		double shaker_freq = argv4; // Set shaker frequency 
//		std::cout << "Vibration frequency:" << std::endl << "Input frequency in Hz: ";
//		std::cin >> shaker_freq;	
		while (!std::cin || shaker_freq <= 0)
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter a value above zero: ";
			std::cin >> shaker_freq;
		}
//		std::cout << std::endl;

		double shaker_amp_temp = argv5; // Set shaker amp
//		std::cout << "Amplitude for vibrations:" << std::endl << "Input amplitude of the particle in mm: ";
//		std::cin >> shaker_amp_temp;
		double shaker_amp = shaker_amp_temp * 1e-3;
		while (!std::cin || shaker_amp <= 0 || shaker_amp > 4*r0)
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please value above 0 that is less that the diameter of 2 particles: ";
			std::cin >> shaker_amp_temp;
			shaker_amp = shaker_amp_temp * 1e-3;
		}
//		std::cout << std::endl;
		
		InfiniteWall basewall; // Initiate infinite wall for the base, lid and x walls
		InfiniteWall lid;
		InfiniteWall xwall;
		PeriodicBoundary ywall; // Initiate a periodic boundary for the y wall, making thw system 2D

		basewall.setSpecies(species);
		lid.setSpecies(species); // Applies the species to the walls
		xwall.setSpecies(species);

		basewall.set(Vec3D( 0.0, 0.0,-1.0),Vec3D(0.0, 0.0, -getZMin())); // Sets initial position of the base
		basewall.setPrescribedPosition([shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating base
		{
			double t = time - 0.10;
			if (t>0.0)
			{
				return Vec3D(0.0, 0.0, shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI));
			}
			else
			{	    
				return Vec3D(0.0,0.0,getZMin());
			}
		});
		wallHandler.copyAndAddObject(basewall); // Adds wall to the system

		lid.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0.0, 0.0, getZMax())); // Sets position of the lid
		wallHandler.copyAndAddObject(lid); // Adds wall to the system
		
		xwall.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(-getXMin(), 0.0, 0.0)); // Sets position of the x min wall
		wallHandler.copyAndAddObject(xwall); // Adds wall to the system

		xwall.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, 0.0)); // Sets  position of the x max wall
		wallHandler.copyAndAddObject(xwall); // Adds wall to the system

		ywall.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax()); // Sets limits for y boundary 
		boundaryHandler.copyAndAddObject(ywall); // Adds boundary to the system

		double w = argv6; // dummy variable for checking user input 
//		std::cout << "Sawtooth Orientation:" << std::endl << "Input 0 for up or 1 for down or 2 for opposite: ";
//		std::cin >> w;
		while (!std::cin || w<0 || w>2 || ceil(w)!=floor(w)) // Conditions to check user input is valid 
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter 0, 1 or 2: ";
			std::cin >> w;
		}

		if (w == 0) // For up orientated sawteeth
		{
			for (float block = 0; block <= 19; block++)
			{
			
				IntersectionOfWalls saw1;
				saw1.setSpecies(species);
		
				saw1.addObject(Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, 0.0, block * 10e-3));
				saw1.addObject(Vec3D(-0.5, 0.0, -0.5), Vec3D(5e-3, 0.0, 5e-3 + (block * 10e-3)));
				saw1.addObject(Vec3D(0.0, 0.0, 1.0), Vec3D(0.0, 0.0, block * 10e-3));

				wallHandler.copyAndAddObject(saw1);

				IntersectionOfWalls saw2;
				saw2.setSpecies(species);
			
				saw2.addObject(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, block * 10e-3));
				saw2.addObject(Vec3D(0.5, 0.0, -0.5), Vec3D(getXMax() - 5e-3, 0.0, 5e-3 + (block * 10e-3)));
				saw2.addObject(Vec3D(0.0, 0.0, 1.0), Vec3D(0.0, 0.0, block * 10e-3));

				wallHandler.copyAndAddObject(saw2);
			}
		}
		else if (w == 1) // For down orientated sawteeth
		{
			for (float block = 0; block <= 19; block++)
			{
			
				IntersectionOfWalls saw1;
				saw1.setSpecies(species);
		
				saw1.addObject(Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, 0.0, block * 10e-3));
				saw1.addObject(Vec3D(-0.5, 0.0, 0.5), Vec3D(5e-3, 0.0, 5e-3 + (block * 10e-3)));
				saw1.addObject(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, (block + 1) * 10e-3));

				wallHandler.copyAndAddObject(saw1);

				IntersectionOfWalls saw2;
				saw2.setSpecies(species);
			
				saw2.addObject(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, block * 10e-3));
				saw2.addObject(Vec3D(0.5, 0.0, 0.5), Vec3D(getXMax() - 5e-3, 0.0, 5e-3 + (block * 10e-3)));
				saw2.addObject(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, (block + 1) * 10e-3));

				wallHandler.copyAndAddObject(saw2);
			}
		}
		else
		{
			for (float block = 0; block <= 19; block++)
			{
			
				IntersectionOfWalls saw1;
				saw1.setSpecies(species);
		
				saw1.addObject(Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, 0.0, block * 10e-3));
				saw1.addObject(Vec3D(-0.5, 0.0, -0.5), Vec3D(5e-3, 0.0, 5e-3 + (block * 10e-3)));
				saw1.addObject(Vec3D(0.0, 0.0, 1.0), Vec3D(0.0, 0.0, block * 10e-3));

				wallHandler.copyAndAddObject(saw1);

				IntersectionOfWalls saw2;
				saw2.setSpecies(species);
			
				saw2.addObject(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, block * 10e-3));
				saw2.addObject(Vec3D(0.5, 0.0, 0.5), Vec3D(getXMax() - 5e-3, 0.0, 5e-3 + (block * 10e-3)));
				saw2.addObject(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, (block + 1) * 10e-3));

				wallHandler.copyAndAddObject(saw2);
			}
		}
		
				int maxpart = (getXMax() - (2.5e-3 * 8)) / (2*r0);

		double area0 = M_PI * r0*r0;
		double area1 = M_PI * r1*r1;
		double totarea = (getXMax() - 20e-3) * getZMax();
		double partarea = 0.0;
//		std::cout << std::endl;
		double breakpercent;
		
//		std::cout << "Density of the system" << std::endl;
//		std::cout << "Very Dilute = 0" << std::endl << "Dilute = 1" << std::endl << "Dense = 2" << std::endl << "Input a value to assign a density: ";
		double d = argv8;		
//		std::cin >> d;
		while (!std::cin || d < 0 || d > 2 || ceil(d)!=floor(d) )
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter 0, 1 or 2: ";
			std::cin >> d;
		}

		if (d == 0)
		{
			breakpercent = 5.0;
		}
		else if (d == 1)
		{
			breakpercent = 10.0;
		}
		else
		{
			breakpercent = 25.0;
		}

		BaseParticle p0;
		p0.setRadius(r0);
		p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
		p0.setSpecies(species);
		
		BaseParticle p1;
		p1.setRadius(r1);
		p1.setVelocity(Vec3D(0.0, 0.0, 0.0));
		p1.setSpecies(species);

		int k = 0;			

		do 
		{ 
			k++;
			
			for (double i = 1.0; i <= maxpart; i++)
			{
				if (B == 1)
				{
					p0.setPosition(Vec3D((i * (2 * r0)) + (7.5e-3), getYMax()/2.0, (k * 2.0 * r0) + (getZMax() * 0.01)));
					particleHandler.copyAndAddObject(p0);
					partarea += area0;
				}

				if (B == 2)
				{
					if (ceil(i/2)==floor(i/2))
					{					
						p0.setPosition(Vec3D((i * (2 * r0)) + (7.5e-3), getYMax()/2.0, (k * 2.0 * r0) + (getZMax() * 0.01)));
						particleHandler.copyAndAddObject(p0);
						partarea += area0;
					}
					else
					{
						p1.setPosition(Vec3D((i * (2 * r0)) + (7.5e-3), getYMax()/2.0, (k * 2.0 * r0) + (getZMax() * 0.01)));
						particleHandler.copyAndAddObject(p1);
						partarea += area1;
					}
				}
			}
		} while ( ((partarea/totarea) * 100) < breakpercent );
		
		BaseParticle wobbler; // Particle to displace main body of particles
		wobbler.setVelocity(Vec3D(-0.5, 0.0, -2.0));
		wobbler.setRadius(r0);
		wobbler.setPosition(Vec3D(getXMax()/2, getYMax()/2, 0.3 * getZMax()));
		wobbler.setSpecies(species); // Makes this bead the same as the rest
		particleHandler.copyAndAddObject(wobbler);
		
	}
};

int main(int argc, char *argv[])
{
	Box Box2DWallStaticTeeth; // Initialise an instance of the above class

	Box2DWallStaticTeeth.setName("Box2DWallStaticTeeth"); // Sets name for output files
	Box2DWallStaticTeeth.setGravity(Vec3D(0.0, 0.0, -9.81)); // Gravity can be set to whatever is needed
	Box2DWallStaticTeeth.setSaveCount(500);

	double tmax = 5.0; // Set run time for simulation
	double tsteps = 5.0e5; // Set number of time steps

//	std::cout << std::endl;
	double b = atof(argv[1]);
//	std::cout << "Single or dual particle system" << std::endl << "Input 1 for single particle or 2 for dual particle system: ";
//	std::cin >> b;
	while (!std::cin || b < 1 || b > 2 || ceil(b)!=floor(b))
	{
		std::cin.clear();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cout << "Please enter 1 or 2: ";
		std::cin >> b;
	}
	Box2DWallStaticTeeth.B = b;

	double r;
	double r_0;
	double r_1;
	double R = 0.0;
	double R0;
	double R1;

	if (b == 1)
	{
//		std::cout << std::endl;	
//		std::cout << "Diameter of particle" << std::endl << "Input diameter of the particle in mm: ";
//		std::cin >> r;
		r = atof(argv[2]);
		R = r * 1e-3 / 2;
		while (!std::cin || R <= 0 || R > 10e-3)
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter a number in mm (maximum 20mm): ";
			std::cin >> r;
			R = r * 1e-3 / 2;
		}
		Box2DWallStaticTeeth.r0 = R;
	}
	else if (b == 2)
	{
//		std::cout << std::endl;
//		std::cout << "Diameter of particle 1, must be smaller than particle 2: " << std::endl << "Input diameter of the particle in mm: ";
//		std::cin >> r_1;
		r_1 = atof(argv[2]);
		R1 = r_1 * 1e-3 / 2;
		while (!std::cin || R1 <= 0 || R1 > 10e-3)
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter a number in mm (maximum 20mm): ";
			std::cin >> r;
			R1 = r_1 * 1e-3 / 2;
		}
			Box2DWallStaticTeeth.r1 = R1;

//		std::cout << std::endl;
//		std::cout << "Diameter of particle 2, must be bigger than particle 1: " << std::endl << "Input diameter of the particle in mm: ";
//		std::cin >> r_0;		
		r_0 = atof(argv[3]);
		R0 = r_0 * 1e-3 / 2;
		while (!std::cin || R0 <= 0 || R0 > 10e-3 || R0 < R1)
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter a number in mm (more than particle 1): ";
			std::cin >> r;
			R0 = r_0 * 1e-3 / 2;
		}
			Box2DWallStaticTeeth.r0 = R0;
	}
	
	if (R == 0.0)
	{
		R = R0;
	}

	Box2DWallStaticTeeth.argv4 = atof(argv[4]);
	Box2DWallStaticTeeth.argv5 = atof(argv[5]);
	Box2DWallStaticTeeth.argv6 = atof(argv[6]);
	Box2DWallStaticTeeth.argv7 = atof(argv[7]);
	Box2DWallStaticTeeth.argv8 = atof(argv[8]);


	Box2DWallStaticTeeth.setXMax(0.12); // Set up the max/min values for the dimensions of the box
	Box2DWallStaticTeeth.setYMax(2 * R);
	Box2DWallStaticTeeth.setZMax(0.2);
	Box2DWallStaticTeeth.setXMin(0.0);
	Box2DWallStaticTeeth.setYMin(0.0);
	Box2DWallStaticTeeth.setZMin(0.0);
	Box2DWallStaticTeeth.setTimeMax(tmax);
	Box2DWallStaticTeeth.setTimeStep(tmax/tsteps);
	
	Box2DWallStaticTeeth.setXBallsAdditionalArguments("-solidf -v0");
    Box2DWallStaticTeeth.solve(argc, argv);
}
