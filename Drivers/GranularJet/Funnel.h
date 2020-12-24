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

#ifndef FUNNEL_H
#define FUNNEL_H
#include "Chute.h"

///Funnel adds a funnel to the Chute class. Parameters needed are:
/// - The number of particles along the sidewall of the funnel
/// - The angle in which the funnel grows
/// - The diameter at the bottom of the funnel
/// - The origin of the funnel: (x,y) coordinates of the middlepoint of the circle
/// - The falling height from the bottom wall
/// - The filling ratio of the funnel, in which (0.0-1.0) range should the funnel be filled

class Funnel : public Chute {
	public:
	
	///This is the default constructor. All it does is set sensible defaults.
	Funnel() {constructor();}
	
	/////Copy-constructor for creates an HGRID problem from an existing MD problem
	//Funnel(MD& other) : MD(other), Mercury3D(other), Chute() {constructor();}
	//Funnel(HGRID_base& other) :  MD(other), Mercury3D(other), Chute() {constructor();}
	//Funnel(Mercury3D& other) : MD(other), Mercury3D(other), Chute() {constructor();}

	///This function prints all funnel data in human readable format
	void write(std::ostream& os, bool writeAllParticles = false) const override;
	
	///This function reads all funnel data
	void read(std::istream& is, ReadOptions opt = ReadOptions::ReadAll) override;
	
	///This is the actual constructor.
	void constructor();
	
	///Here we define the outflow
	void cleanChute();
	
	///Set funnel number of particles along the funnel heigth:
	void set_funnz(double funnz_){funnz=funnz_;}
	
	///Get the funnel number of particles along the funnel heigth:
	double get_funnz() const {return funnz;}
	
	///Set the filling ratio
	void set_funfr(double funfr_){
	  if (funfr_<0.0 || funfr_>1.0) {std::cerr << "Filling Ratio is below 0 or above 1, default value of 0.33 is used!"<<std::endl; }
		else {funfr=funfr_;}
	}
	
	///Get the filling ratio
	double get_funfr() const {return funfr;}
	
	///Get funnel radius:
	double get_funr(){return funr;}
	
	///Set funnel origin of the funnel:
	void set_funO(double x, double y){funO[0]=x; funO[1]=y;}
	void set_funO(double* x){funO[0]=x[0]; funO[1]=x[1];}
	
	///Get funnel origin of the funnel:
	const double* get_funO() const {return funO;}
	
	double get_funOx() const {return funO[0];}
	
	double get_funOy() const {return funO[1];}
	
	double get_fundiag() const {return fundiag;}
	
	double get_funrmax() const {return funrmax;}
	
	///Set funnel angle:
	void set_funa(double funa_){funa=funa_*constants::pi/180;}
	
	///Get funnel angle:
	double get_funa() const {return funa;}
		
	///Get funnel Heigth:
	double get_funH() const {return funH;}
	
	///Set falling heigth:
	void set_funHf(double funHf_){funHf=funHf_;}
	
	///Get falling heigth:
	double get_funHf() const {return funHf;}
	
	///Set minimum funnel diameter:
	void set_funD(double funD_){funD=funD_;}
	
	///Get minimum funnel diameter:
	double get_funD() const {return funD;}
	
	void setName_();

protected:
	SphericalParticle inflowParticle_;
	
	///Set funnel Heigth:
	void set_funH(double funH_){funH=funH_;}

	void set_funrmax(double funrmax_){funrmax=funrmax_;}

	void set_fundiag(double fundiag_){fundiag=fundiag_;}

	///Set funnel radius:
	void set_funr(double funr_){funr=funr_;}

	/// initialise particle position, velocity, radius
	virtual void setupInitialConditions() override;
	
	///Sets variable values for particles that are created at the inflow
	virtual void create_inflow_particle();
	
	///Create the funnel
	virtual void create_funnel();
	
	///Updates the parameters for the funnel;
	virtual void update_funnel();
	
	///Check the funnel parameters
	virtual void check_funnel();
	
	///Create or update the walls
	virtual void create_walls();
	
	virtual bool readNextArgument(int& i, int argc, char *argv[]) override;

	double funr; // Funnel radius.
	double funO[2]; // Origin of the funnel is the (x,y) location of the center of the top of the funnel.
	double funa; // Angle of the funnel
	double funH; // Heigth of the funnel
	double funHf; // Falling heigth
	double funD; // Funnel diameter at the downside of the funnel.
	double funnz; //Number of particles along the heigth of the funnel
	double funfr; // Filling ratio of the funnel
	double fundiag; //The diagonal of the filling region
	double funrmax; //The maximum range for r
	
};
#endif
