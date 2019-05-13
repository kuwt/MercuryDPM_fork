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
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include <math.h>

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
		auto species = new LinearViscoelasticSpecies;
	
		speciesHandler.addObject(species);
		setParticleDimensions(3);
	
		double YsMod = 3.0; // Young Modulus of Nylon
		//Supposed to be 3 GPa = 3e9 Pa??

//		std::cout << std::endl;

		double stiffness = (M_PI/2.0) * r0 * YsMod * 100 * 50; // Stiffness of Nylon
		// Need proper equation for this

		double density = 1150.0; // Denisty of Nylon
		double mass = M_PI * density * r0*r0*r0 * 4.0 / 3.0;

		species -> setDensity(density);
		species -> setStiffnessAndRestitutionCoefficient(stiffness, 0.9, mass);
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
// 		std::cout << std::endl;

		double shaker_amp_temp = argv5; // Set shaker amp
//		std::cout << "Amplitude fo vibrations:" << std::endl << "Input amplitude of the particle in mm: ";
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
		InfiniteWall sidewall;

		basewall.setSpecies(species);
		lid.setSpecies(species);
		sidewall.setSpecies(species);

		basewall.set(Vec3D( 0.0, 0.0, -1.0),Vec3D(0.0, 0.0, -getZMin()));
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
		wallHandler.copyAndAddObject(basewall);

		lid.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0.0, 0.0, getZMax()));
		wallHandler.copyAndAddObject(lid);
		
		sidewall.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(-getXMin(), 0.0, 0.0));
		wallHandler.copyAndAddObject(sidewall);

		sidewall.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, 0.0));
		wallHandler.copyAndAddObject(sidewall);

		sidewall.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0, -getYMin(), 0.0));
		wallHandler.copyAndAddObject(sidewall);

		sidewall.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0, getYMax(), 0.0));
		wallHandler.copyAndAddObject(sidewall);

		SphericalParticle saw; // Fixed particles to make up the sawteeth
		saw.setRadius(1.25e-3);
		saw.setVelocity(Vec3D(0.0, 0.0, 0.0));
		saw.setSpecies(species);
		
		double w = argv6;
//		std::cout << "Sawtooth Orientation:" << std::endl << "Input 0 for down or 1 for up or 2 for opposite corners: ";
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
			std::cout << "Please enter 0, 1: ";
			std::cin >> v;
		}

		for (float block = 0; block <= 19; block++)
		{
			for (float k = 0; k <= 3; k++) for (float j = 0; j <= 47; j++) for (float i = 0; i <= 47; i++)
			{
				if (w == 0) // Down orienatated sawteeth			
				{
					if (k >= i && i <= 3) // Sawteeth on X walls, down orientation
					{
						saw.setPosition(Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, j * 2.5e-3, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					
						saw.setPosition(Vec3D(getXMax() - (i * 2.5e-3), j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(getXMax() - (i * 2.5e-3), j * 2.5e-3, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(getXMax() - (i * 2.5e-3), j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}

					if (k >= j && j <= 3) // Sawteeth on Y walls, down orientation
					{
						saw.setPosition(Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, j * 2.5e-3, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					
						saw.setPosition(Vec3D(i * 2.5e-3, getYMax() - (j * 2.5e-3), (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, getYMax() - (j * 2.5e-3), ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(i * 2.5e-3, getYMax() - (j * 2.5e-3), (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}
				}
				if (w == 1) // Up orienatated sawteeth			
				{
					if (k + i <= 3 && i <= 3) // Sawteeth on X walls, up orientation
					{
						saw.setPosition(Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, j * 2.5e-3, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
				
						saw.setPosition(Vec3D(getXMax() - (i * 2.5e-3), j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(getXMax() - (i * 2.5e-3), j * 2.5e-3, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(getXMax() - (i * 2.5e-3), j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}
					if (k + j <= 3 && j <= 3) // Sawteeth on Y walls, up orientation
					{
						saw.setPosition(Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 + (k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, j * 2.5e-3, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 + (k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
				
						saw.setPosition(Vec3D(i * 2.5e-3, getYMax() - (j * 2.5e-3), (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, getYMax() - (j * 2.5e-3), ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(i * 2.5e-3, getYMax() - (j * 2.5e-3), (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}
				}
					
				else // Opposite orienatated sawteeth			
				{
					if (k <= i && i <= 3) // Sawteeth on X walls, down orientation
					{
						saw.setPosition(Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, j * 2.5e-3, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
				
						saw.setPosition(Vec3D(getXMax() - (i * 2.5e-3), j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, k, j, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(getXMax() - (i * 2.5e-3), j * 2.5e-3, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(getXMax() - (i * 2.5e-3), j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}

					if (k + j <= 3 && j <= 3) // Sawteeth on Y walls, up orientation
					{
						saw.setPosition(Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, j * 2.5e-3, ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(i * 2.5e-3, j * 2.5e-3, (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
				
						saw.setPosition(Vec3D(i * 2.5e-3, getYMax() - (j * 2.5e-3), (1.25 +(k * 2.5) + (block * 10))*1e-3)); 
						saw.fixParticle();
						if (v == 0)
						{						
							saw.setPrescribedPosition([i, j, k, block, shaker_freq, shaker_amp, this] (double time) // Function to cause vibrating sawteeth
							{
								double t = time - 0.05;
								if (t > 0.0)
								{
									return Vec3D(i * 2.5e-3, getYMax() - (j * 2.5e-3), ((1.25 + (k * 2.5)+(block * 10)) * 1e-3)+(shaker_amp * std::sin(t * 2.0 * shaker_freq * M_PI)));
								}
								else
								{	    
									return Vec3D(i * 2.5e-3, getYMax() - (j * 2.5e-3), (1.25 +(k * 2.5) + (block * 10))*1e-3);
								}
							});
						}
						particleHandler.copyAndAddObject(saw);
					}
				}
			}
		}

		int maxpart = (getXMax() - (2.5e-3 * 8)) / (2*r0);

		double vol0 = 4/3 * M_PI * r0*r0*r0;
		double vol1 = 4/3 * M_PI * r1*r1*r1;
		double totvol = (getXMax() - 20e-3) * (getYMax() - 20e-3) * getZMax();
		double partvol = 0.0;
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
			breakpercent = 1.5;
		}
		else if (d == 1)
		{
			breakpercent = 4.0;
		}
		else
		{
			breakpercent = 10.0;
		}

		SphericalParticle p0;
		p0.setRadius(r0);
		p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
		p0.setSpecies(species);
		
		SphericalParticle p1;
		p1.setRadius(r1);
		p1.setVelocity(Vec3D(0.0, 0.0, 0.0));
		p1.setSpecies(species);

		int k = 0;			

		do 
		{ 
			k++;
			
			for (double i = 1.0; i <= maxpart; i++) for (double j = 1.0; j <= maxpart; j++)
			{
				if (B == 1)
				{
					p0.setPosition(Vec3D((i * (2 * r0)) + (7.5e-3), (j * (2 * r0)) + (7.5e-3), (k * 2.0 * r0) + (getZMax() * 0.01)));
					particleHandler.copyAndAddObject(p0);
					partvol += vol0;
				}

				if (B == 2)
				{
					if (ceil(i/2)==floor(i/2))
					{					
						p0.setPosition(Vec3D((i * (2 * r0)) + (7.5e-3), (j * (2 * r0)) + (7.5e-3), (k * 2.0 * r0) + (getZMax() * 0.01)));
						particleHandler.copyAndAddObject(p0);
						partvol += vol0;
					}
					else
					{
						p1.setPosition(Vec3D((i * (2 * r0)) + (7.5e-3), (j * (2 * r0)) + (7.5e-3), (k * 2.0 * r0) + (getZMax() * 0.01)));
						particleHandler.copyAndAddObject(p1);
						partvol += vol1;
					}
				}
			}
		} while (((partvol/totvol) * 100) < breakpercent);
		
		SphericalParticle wobbler;
		wobbler.setVelocity(Vec3D(-0.1, 0.0, -1.0));
		wobbler.setRadius(r0);
		wobbler.setPosition(Vec3D(0.55 * getXMax(), getYMax()/2, 0.75 * getZMax()));
		wobbler.setSpecies(species);
		particleHandler.copyAndAddObject(wobbler);
		
	}
};

int main(int argc, char *argv[])
{
	Box Box3DPartTeeth;

	Box3DPartTeeth.setName("Box3DPartTeeth");
	Box3DPartTeeth.setGravity(Vec3D(0.0, 0.0, -9.81));
	Box3DPartTeeth.setSaveCount(500);

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
	Box3DPartTeeth.B = b;

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
		Box3DPartTeeth.r0 = R;
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
			Box3DPartTeeth.r1 = R1;

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
			Box3DPartTeeth.r0 = R0;
	}
	
	if (R == 0.0)
	{
		R = R0;
	}

	Box3DPartTeeth.argv4 = atof(argv[4]);
	Box3DPartTeeth.argv5 = atof(argv[5]);
	Box3DPartTeeth.argv6 = atof(argv[6]);
	Box3DPartTeeth.argv7 = atof(argv[7]);
	Box3DPartTeeth.argv8 = atof(argv[8]);

	Box3DPartTeeth.setXMax(0.12);
	Box3DPartTeeth.setYMax(0.12);
	Box3DPartTeeth.setZMax(0.2);
	Box3DPartTeeth.setXMin(0.0);
	Box3DPartTeeth.setYMin(0.0);
	Box3DPartTeeth.setZMin(0.0);
	Box3DPartTeeth.setTimeMax(tmax);
	Box3DPartTeeth.setTimeStep(tmax/tsteps);
	
	Box3DPartTeeth.setXBallsAdditionalArguments("-solidf -v0");

    Box3DPartTeeth.solve(argc, argv);
}
