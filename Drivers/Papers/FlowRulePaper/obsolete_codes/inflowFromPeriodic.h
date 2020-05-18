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

#include "scr/Chute.h"

class inflowFromPeriodic : public Chute
{
public:
	inflowFromPeriodic() {
		//load particles and chute setup into periodic chute
		setName("../../H20A22L0.5M0.5");
		load_restart_data();
		setRestarted(false);

		// adds new species
		speciesHandler.copyAndAddObject(speciesHandler.getObject(0));

		// creates bottom outside periodic chute of species 1
		set_Nmax(get_N()*12.);
		int N=get_N();
		for (int j=0; j<5; j++) {
			for (int i=0; i<N; i++) {
				if (Particles[i].is_fixed()) {
					Particles.push_back(Particles[i]);
					Particles.back().indSpecies=1;
					Particles.back().Position.X+=20*j;
				}
			}
		}
		
		Chute SlowBottom;
		SlowBottom.setName("../../H20A22L0.75M0.5");
		SlowBottom.load_restart_data();
		for (int j=5; j<10; j++) {
			for (int i=0; i<SlowBottom.get_N(); i++) {
				if (SlowBottom.getObjects()[i].is_fixed()) {
					Particles.push_back(SlowBottom.getObjects()[i]);
					Particles.back().indSpecies=1;
					Particles.back().Position.X+=20*j;
				}
			}
		}
		setXMax(10*getXMax());

		Walls.resize(0);
	}

	inflowFromPeriodic(string restart_file) {
		//load particles and chute setup into periodic chute
		setName(restart_file.c_str());
		load_restart_data();
		setName("restarted");
		set_HGRID_num_buckets_to_power();
	}

	//Do not add or remove particles
	void actionsBeforeTimeStep(){cleanChute();};

	void cleanChute() {
		//clean outflow every 100 time steps
		static int count = 0, maxcount = 100;
		if (count>maxcount)
		{
			count = 0;
			// delete all outflowing particles
			for (unsigned int i=0;i<Particles.size();)
			{
				if (Particles[i].Position.Z<getZMin()-10)//||Particles[i].Position.Z+Particles[i].Radius<zmin)
				{
					#ifdef DEBUG_OUTPUT_FULL
						cout << "erased:" << Particles[i] << endl;
					#endif
					removeParticle(i);
				}	
				else i++;
			}
		} else count++;
	}
	
	//Do not add bottom
	void setupInitialConditions(){}
	
	void integrateBeforeForceComputation(int i)
	{
		Particles[i].Velocity += Particles[i].Force*Particles[i].get_invmass()*0.5*getTimeStep();
		Particles[i].Position += Particles[i].Velocity*getTimeStep();
		if (speciesHandler.getObject(0)->mu)
		{	
			Particles[i].AngularVelocity += Particles[i].Torque*Particles[i].get_invinertia()*0.5*getTimeStep();
			Particles[i].Angle += Particles[i].AngularVelocity*getTimeStep();
		}
		
		// This shifts particles that moved through periodic walls
		if (Particles[i].indSpecies==0 && WallsPeriodic[0].distance(Particles[i])<0) {
			if (Particles[i].Position.X>18 && Particles[i].Position.X<22) {
				if (get_Nmax()<=get_N()) set_Nmax(get_Nmax()+10000);
				Particles.push_back(Particles[i]);
				Particles.back().indSpecies=1;
			}
			WallsPeriodic[0].shift_position(Particles[i].Position);
		}
		if (WallsPeriodic[1].distance(Particles[i])<0)
			WallsPeriodic[1].shift_position(Particles[i].Position);
	}

	//~ void computeExternalForces(int CI) {
		//~ MD::computeExternalForces(CI);
		//limit force
		//~ if (Particles[CI].indSpecies==1 && Particles[CI].Position.X<40) {
			//~ Particles[CI].Force *= max(0.,(Particles[CI].Position.X-22.)/18.);
			//~ if (std::max(0.,(Particles[CI].Position.X-22.)/18.)) 
				//~ cout << max(0.,(Particles[CI].Position.X-22.)/18.) << endl;
		//~ }
	//~ }

	void printTime() {
		int counter=0;
		for (int i=0; i<get_N(); i++)
			if (Particles[i].indSpecies==0) counter++;

		cout << "t=" << setprecision(3) << left << setw(6) << getTime() 
			<< ", Nmax=" << setprecision(3) << left << setw(6) << get_Nmax()
			<< ", counter=" << setprecision(3) << left << setw(6) << counter << endl;
	}

	void add_forces(int PI, int PJreal) {
			// Add and subtract this force to the force acting on particles i, j
			if (Particles[PI    ].indSpecies!=Particles[PJreal].indSpecies) {
				if (Particles[PI    ].indSpecies==0) {
					if (Particles[PI    ].Position.X>18) {
						Particles[PJreal].Force-=force;
					}
				} else {
					if (Particles[PJreal].Position.X>18) {
						Particles[PI    ].Force+=force;					
					}
				}
			} else {
				Particles[PI    ].Force+=force;
				Particles[PJreal].Force-=force;				
			}		
	}

	int Check_and_Duplicate_Periodic_Particle(int i, int nWallPeriodic)
	{
		int C=0; //Number of particles created
		if (nWallPeriodic==0 && Particles[i].indSpecies==0 && WallsPeriodic[0].distance(Particles[i])<Particles[i].Radius+get_radius_of_largest_particle())
		{
			CParticle F0=Particles[i];
			WallsPeriodic[0].shift_position(F0.Position);
			
			//If Particle is double shifted, get correct original particle
			int From=i;
			while(Particles[From].get_periodicFromParticle()!=-1)
				From=Particles[From].get_periodicFromParticle();		
			F0.set_periodicFromParticle(From);
									
			Particles.push_back(F0);
			HGridInsertParticleToHgrid(Particles[get_N()-1]);
			C++;
			
			//Check for double shifted particles
			C+=Check_and_Duplicate_Periodic_Particle(get_N()-1, 0+1);
		}
		for (int k=max(1,nWallPeriodic); k<get_NWallPeriodic(); k++) {	//Loop over all still posible walls
			///\todo{distance should be Particles[i].Radius+max(Particles[j].Radius)}
			if (WallsPeriodic[k].distance(Particles[i])<Particles[i].Radius+get_radius_of_largest_particle())
			{
				CParticle F0=Particles[i];
				WallsPeriodic[k].shift_position(F0.Position);
				
				//If Particle is double shifted, get correct original particle
				int From=i;
				while(Particles[From].get_periodicFromParticle()!=-1)
					From=Particles[From].get_periodicFromParticle();		
				F0.set_periodicFromParticle(From);
										
				Particles.push_back(F0);
				HGridInsertParticleToHgrid(Particles[get_N()-1]);
				C++;
				
				//Check for double shifted particles
				C+=Check_and_Duplicate_Periodic_Particle(get_N()-1, k+1);
			}
		}
		return(C);
	}

	void writeXBallsScript()
	{
		
		stringstream file_name;
		ofstream script_file;
		file_name << problem_name.str() <<".disp";
		script_file.open((file_name.str()).c_str());

		///First put in all the script lines. All these lines do is move you to the correct directory from any location
		script_file << "#!/bin/bash" << endl;
		script_file << "x=$(echo $0 | cut -c2-)" << endl;
		script_file << "file=$PWD$x" << endl;
		script_file << "dirname=`dirname \"$file\"`" << endl;
		script_file << "cd $dirname" << endl;
		
		double scale;
		int format;

		int width=1700-140, height=1000-140;
		int ratio=getSystemDimensions()<3?(getXMax()-getXMin())/(getYMax()-getYMin()):(getXMax()-getXMin())/(getZMax()-getZMin());
		if  (ratio>width/height) height=width/ratio;
		else width=height*ratio;

		if (getSystemDimensions()<3) 
			{ // dim = 1 or 2
			format = 8;
			if (getXBallsScale()<0)
				{
				scale = 1.0 / max( getYMax()-getYMin(), getXMax()-getXMin() );
				}
			else
				{
				scale=getXBallsScale();
				}
			} 
		else 
			{ //dim==3
			format = 14;
			if (getXBallsScale()<0)
				{
					scale = width/480/(getXMax()-getXMin());
				}
	
			else
				{
					scale=getXBallsScale();
				}
			
			}
					
		script_file << "../xballs -format " << format 
			<< " -f " << data_filename.str() << ((get_options_data()==2)?"'.0000":"")
			<< " -s " << scale
			<< " -w " << width+140
			<< " -h " << height+140
			<< " -cmode " << getXBallsColourMode()
			<< " -cmax -scala 4 -sort " 
			<< getXBallsAdditionalArguments()
			<< " $*";
		if (getXBallsVectorScale()>-1)
			{
				script_file << " -vscale " << getXBallsVectorScale();
			}
		script_file.close();
		
		//This line changes teh file permision and give the owener (i.e. you) read, write and excute permission to the file.
		chmod((file_name.str().c_str()),S_IRWXU);

	}


	//~ void outputXBallsData()
	//~ {
		//~ format=14;
		//~ // This outputs the location of walls and how many particles there are to file this is required by the xballs plotting
		//~ int counter=0;
		//~ for (unsigned int i = 0;i<Particles.size();i++) if (Particles[i].indSpecies==0) counter++;
		//~ if (format!=14) // dim = 1 or 2
			//~ data_file << Particles.size()-counter << " " <<t << " " 
			//~ << xmin << " " << ymin << " " 
			//~ << xmax << " " << ymax << " " << endl;
		//~ else //dim==3
			//~ data_file << Particles.size()-counter << " " <<t << " " 
			//~ << xmin << " " << ymin << " " << zmin << " "
			//~ << xmax << " " << ymax << " " << zmax << " " << endl;
					//~ // This outputs the particle data
		//~ for (unsigned int i = 0;i<Particles.size();i++) {
			//~ if (Particles[i].indSpecies==1) outputXBallsDataParticle(i);
		//~ }
	//~ }

};
