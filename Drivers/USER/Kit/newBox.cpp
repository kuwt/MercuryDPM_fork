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
#include "Boundaries/PeriodicBoundary.h"
#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
class Box : public Mercury3D
{
public:
	
	double r0;
	double r1;
	int B;
	double argv4;
	double argv5;
	int argv6;
	double argv7;
	double argv8;

	double toothRestitution;

	void setupInitialConditions()
	{
		//general particle species
		auto species = new LinearViscoelasticSpecies;	
		speciesHandler.addObject(species);
		//sawtooth particle species
		auto toothSpecies = new LinearViscoelasticSpecies; 
		speciesHandler.addObject(toothSpecies);
		setParticleDimensions(3);
	
		double YsMod = 3.0; // Young Modulus of Nylon
		//Supposed to be 3 GPa = 3e9 Pa??
		
//		std::cout << std::endl;

		double stiffness = (M_PI/2.0) * r0 * YsMod * 100 * 50; // Stiffness of Nylon
		// Need proper equation for this

		double density = 1150.0; // Denisty of Nylon
		double mass = M_PI * density * r0*r0*r0 * 4.0 / 3.0;

		//setting particle species properties
		species -> setDensity(density);
		species -> setStiffnessAndRestitutionCoefficient(stiffness, 0.9, mass);
		
		//setting tooth particle species properties
		toothSpecies -> setDensity(density);
		toothSpecies -> setStiffnessAndRestitutionCoefficient(stiffness, toothRestitution, mass);

		//setting mixed interaction properties
		speciesHandler.getMixedObject(species, toothSpecies)->setStiffnessAndRestitutionCoefficient(stiffness, toothRestitution, 2*mass);
		
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
		
		InfiniteWall basewall;
		InfiniteWall lid;
		InfiniteWall xwall;
		PeriodicBoundary ywall;		

		basewall.setSpecies(species);
		lid.setSpecies(species);
		xwall.setSpecies(species);

		basewall.set(Vec3D( 0.0, 0.0,-1.0),Vec3D(0.0, 0.0, -getZMin()));
		basewall.setPrescribedPosition([shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating base
		{
			double t = time - 0.05;
			if (t>0.0)
			{
				return Vec3D(0.0, 0.0, shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI));
			}
			else
			{	    
				return Vec3D(0.0,0.0,getZMin());
			}
		});
		wallHandler.copyAndAddObject(basewall);

		lid.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0.0, 0.0, getZMax()));
		wallHandler.copyAndAddObject(lid);
		
		xwall.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(-getXMin(), 0.0, 0.0));
		wallHandler.copyAndAddObject(xwall);

		xwall.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, 0.0));
		wallHandler.copyAndAddObject(xwall);

		ywall.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
		boundaryHandler.copyAndAddObject(ywall);

		BaseParticle saw; // Fixed particles to make up the sawteeth
		saw.setRadius(1.25e-3);
		saw.setVelocity(Vec3D(0.0, 0.0, 0.0));
		saw.setSpecies(toothSpecies);

		double w = argv6;
//		std::cout << "Sawtooth Orientation:" << std::endl << "Input 0 for down or 1 for up or 2 for opposite: ";
//		std::cin >> w;
		while (!std::cin || w<0 || w>2 || ceil(w)!=floor(w) )
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter 0, 1 or 2: ";
			std::cin >> w;
		}
//		std::cout << std::endl;

		double v = argv7;
//		std::cout << "Vibrating Sawteeth:" << std::endl << "Input 0 for yes or 1 for no: ";
//		std::cin >> v;
		while (!std::cin || v<0 || v>1 || ceil(v)!=floor(v) )
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter 0 or 1: ";
			std::cin >> v;
		}

		for (float block = 0; block <= 19; block++) // Repeating sawtooth block
		{
			for (float k = 0; k <= 3; k++) for (float i = 0; i <= 3; i++) // Placing each particle in sawtooth
			{
				if (w == 0) // Down orientated sawteeth
				{				
					if (k >= i) // Condition for down orientation
					{
						saw.setPosition(Vec3D(i * 2.5e-3, getYMax()/2, (1.25 + (k * 2.5) + (block * 10)) * 1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, getYMax()/2, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(Vec3D(i * 2.5e-3, getYMax()/2, (1.25 +(k * 2.5) + (block * 10))*1e-3));
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					
						saw.setPosition(Vec3D(getXMax() - (i * 2.5e-3), getYMax()/2, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{
							saw.setPrescribedPosition([i, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(getXMax() - (i * 2.5e-3), getYMax()/2, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(Vec3D(getXMax() - (i * 2.5e-3), getYMax()/2, (1.25 +(k * 2.5) + (block * 10)) * 1e-3));
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}
				}
				else if (w == 1) // Up orientated sawteeth
				{
					if (k + i <= 3) // Condition for up orientation
					{
						saw.setPosition(Vec3D(i * 2.5e-3, getYMax()/2, (1.25 + (k * 2.5) + (block * 10)) * 1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{
							saw.setPrescribedPosition([i, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, getYMax()/2, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(Vec3D(i * 2.5e-3, getYMax()/2, (1.25 +(k * 2.5) + (block * 10))*1e-3));
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					
						saw.setPosition(Vec3D(getXMax() - (i * 2.5e-3), getYMax()/2, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{
							saw.setPrescribedPosition([i, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(getXMax() - (i * 2.5e-3), getYMax()/2, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(Vec3D(getXMax() - (i * 2.5e-3), getYMax()/2, (1.25 +(k * 2.5) + (block * 10)) * 1e-3));
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}
				}
				else
				{
					if (k >= i) // Condition for down orientation
					{
						saw.setPosition(Vec3D(i * 2.5e-3, getYMax()/2, (1.25 + (k * 2.5) + (block * 10)) * 1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{
							saw.setPrescribedPosition([i, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{	
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, getYMax()/2, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(Vec3D(i * 2.5e-3, getYMax()/2, (1.25 +(k * 2.5) + (block * 10))*1e-3));
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}
					if (k + i <= 3) // Condition for up orientation
					{
						saw.setPosition(Vec3D(getXMax() - (i * 2.5e-3), getYMax()/2, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{
							saw.setPrescribedPosition([i, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(getXMax() - (i * 2.5e-3), getYMax()/2, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(Vec3D(getXMax() - (i * 2.5e-3), getYMax()/2, (1.25 +(k * 2.5) + (block * 10)) * 1e-3));
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}
				}
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
		breakpercent = d;

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
		
		BaseParticle wobbler;
		wobbler.setVelocity(Vec3D(-0.1, 0.0, -1.0));
		wobbler.setRadius(r0);
		wobbler.setPosition(Vec3D(0.55 * getXMax(), getYMax()/2, 0.75 * getZMax()));
		wobbler.setSpecies(species);
		particleHandler.copyAndAddObject(wobbler);	
	}
};

int main()
{
	Box Box2DPartTeeth;

	//parameters
	double radS = 3.0; //small particle radius
	double radL = 5.0; //large particle radius
	int poly = 1; //the degree of polydispersity (1 for mono, 2 for bi)
	double freq = 75; //vibrational frequency in Hz
	double amp = 1.104; //vibrational amplitude in mm;
	int orientation = 1; //orientation of sawteeth - 0 for down, 1 for up, 2 for opposing
	bool teethOff = false; //if true, floor vibrates but teeth do not
	double packFrac = 10.0; //the fraction of the box volume occupied by particles; 10.0 gives realistic system
	double eTeeth = 0.9; //the restitution coefficient of the wall particles

	//creating a suitably informative name
	std::stringstream ss;
	std::string baseName = "2teeth";
	ss << baseName << "_" << poly << "species_" << "orientation" << orientation << "_teethOff" << teethOff << "_eWall" << eTeeth; 


	Box2DPartTeeth.setName(ss.str());
	Box2DPartTeeth.setGravity(Vec3D(0.0, 0.0, -9.81));
	Box2DPartTeeth.setSaveCount(1000);

	double tmax = 500.0; // Set run time for simulation
	double timeStep = 1.0 / (50.0 * 800.0); // Set number of time steps

	Box2DPartTeeth.B = poly;

	double r;
	double r_0;
	double r_1;
	double R = 0.0;
	double R0;
	double R1;

	if (poly == 1)
	{
		r = radS;
		R = r * 1e-3 / 2;
		while (!std::cin || R <= 0 || R > 10e-3)
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter a number in mm (maximum 20mm): ";
			std::cin >> r;
			R = r * 1.0e-3 / 2.0;
		}
		Box2DPartTeeth.r0 = R;
	}
	else if (poly == 2)
	{
		r_1 = radS;
		R1 = r_1 * 1e-3 / 2;
		while (!std::cin || R1 <= 0 || R1 > 10e-3)
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter a number in mm (maximum 20mm): ";
			std::cin >> r;
			R1 = r_1 * 1e-3 / 2;
		}
			Box2DPartTeeth.r1 = R1;

//		std::cout << std::endl;
//		std::cout << "Diameter of particle 2, must be bigger than particle 1: " << std::endl << "Input diameter of the particle in mm: ";
//		std::cin >> r_0;		
		r_0 = radL;
		R0 = r_0 * 1e-3 / 2;
		while (!std::cin || R0 <= 0 || R0 > 10e-3 || R0 < R1)
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Please enter a number in mm (more than particle 1): ";
			std::cin >> r;
			R0 = r_0 * 1e-3 / 2;
		}
			Box2DPartTeeth.r0 = R0;
	}
	
	if (R == 0.0)
	{
		R = R0;
	}

	Box2DPartTeeth.argv4 = freq;
	Box2DPartTeeth.argv5 = amp;
	Box2DPartTeeth.argv6 = orientation; 
	Box2DPartTeeth.argv7 = teethOff;
	Box2DPartTeeth.argv8 = packFrac;

	Box2DPartTeeth.setXMax(0.12);
	Box2DPartTeeth.setYMax(2 * R);
	Box2DPartTeeth.setZMax(0.2);
	Box2DPartTeeth.setXMin(0.0);
	Box2DPartTeeth.setYMin(0.0);
	Box2DPartTeeth.setZMin(0.0);
	Box2DPartTeeth.setTimeMax(tmax);
	Box2DPartTeeth.setTimeStep(timeStep);
	
	Box2DPartTeeth.toothRestitution = eTeeth;
	
	Box2DPartTeeth.setXBallsAdditionalArguments("-solidf -v0");

  	Box2DPartTeeth.solve();
}
