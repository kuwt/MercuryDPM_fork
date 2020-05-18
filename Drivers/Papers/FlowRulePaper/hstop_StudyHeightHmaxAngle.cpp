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

/** The files silbert_hstop.cpp and silbert_hstop_isAbove.cpp together define a simple bisection algorithm to determine the hstop curve, i.e. at which angle and height a steady chute flow is not possible anymore and the flow stops.
 * silbert_hstop.cpp starts at a low height and determines the highest angle at which the chute flow subsides, then increases height and repeats the process (it takes into account that the hstop(theta) curve is monotonically decreasing). It does this by calling silbert_hstop_isAbove.cpp, which runs a simulation of a fixed height and angle and returns 0 is the flow subsides before tmax and 1 if not
 **/

#include "SilbertPeriodic.h"

class SilbertHstop : public SilbertPeriodic 
{
public:

	SilbertHstop() {
		pointIsAboveCurve = true;
		//~ setTimeStep(getTimeStep()/3.);
		setSaveCount(1/getTimeStep()); 
		setTimeMax(500);
		fStatFile.setFileType(FileType::NO_FILE);//fStatFile.setFileType(FileType::NO_FILE);
		dataFile.setFileType(FileType::NO_FILE);//dataFile.setFileType(FileType::NO_FILE);
	}

	double pointIsAboveCurve;

	void actionAfterTimeStep()
	{
		SilbertHstop::actionAfterTimeStep();
		double kineticEnergy=getKineticEnergy();
		std::cout << "ene_kin/ene_ele=" << kineticEnergy/getElasticEnergy() << ", piac=" << pointIsAboveCurve << std::endl;
		if (kineticEnergy/getElasticEnergy()<1e-5 && getTime()>10) pointIsAboveCurve=false;
	}

	bool continueSolve() const override {
		if (ceil(getTime())!=ceil(getTime()+getTimeStep())) printTime();
		if (pointIsAboveCurve) return true;
		else return false;
	}

	bool IsAboveCurve() {
		return IsAboveCurve(get_H(),getChuteAngle()/constants::pi*180.);
	}

	bool IsAboveCurve(Mdouble h, Mdouble a) {
		Mdouble func = 40*(tan(a/180.*constants::pi)-tan(30./180.*constants::pi))/(tan(20./180.*constants::pi)-tan(30./180.*constants::pi));
		//~ cout << " h=" << h << ", a=" << a << ", f=" << (h>func) << endl;
		//~ cerr << " h=" << h << ", a=" << a << ", f=" << (h>func) << endl;
		if (h>func) return true;
		else return false;	
	}

	void solve_analytic() {
		setRoughBottomType(MONOLAYER_DISORDERED);
		setTimeMax(1e-4);
		solve();
		pointIsAboveCurve=(tan(getChuteAngle())<.2+.2*exp(getZMax()));
		pointIsAboveCurve=(getZMax()>80-4*getChuteAngleDegrees());
	};

};

std::fstream ReportFile;

bool PointIsAboveCurve(int argc, char *argv[], Mdouble h, Mdouble a, int study_num) {
  std::cout << "PointIsAboveCurve(" << study_num << ", h=" << h << ", a=" << a << ")" << std::endl;
	SilbertHstop problem;
	problem.setInflowHeight(h) ;
	problem.setChuteAngle(a);
	problem.set_study(study_num);
	std::stringstream name;
	name << "H" << problem.getInflowHeight() 
		 << "A" << problem.getChuteAngleDegrees() 
		 << "L" << round(100.*problem.getFixedParticleRadius()*2.)/100.
	     << "M" << problem.species->getSlidingFrictionCoefficient()
		 << "B" << problem.getSlidingFrictionCoefficientBottom();
	problem.setName(name.str().c_str());
	std::cout << "starting " << name.str().c_str() << std::endl;
	problem.writeRestartFile();
	if ((argc>5) && (!strcmp(argv[5],"-test"))) {
		// for test runs:
		problem.solve_analytic();
	} else {
		// for real runs:
		problem.solve();
	}
	std::stringstream com;
	ReportFile   << problem.getInflowHeight() 
		 << "\t" << problem.getChuteAngleDegrees() 
		 << "\t" << round(100.*problem.getFixedParticleRadius()*2.)/100.
		 << "\t" << problem.species->getSlidingFrictionCoefficient()
		 << "\t" << problem.getSlidingFrictionCoefficientBottom()
		 << "\t" << problem.pointIsAboveCurve
		     << std::endl;
	return problem.pointIsAboveCurve;
}

int HstopCurve(int argc, char *argv[], int study_num, Mdouble h, Mdouble hMax, Mdouble a)
{
  std::cout << "restart at study_num=" << study_num << ", h=" << h << ", a=" << a << std::endl;
	
	std::stringstream name;
	name << "Report" << study_num;
	ReportFile.open(name.str().c_str(),std::fstream::out | std::fstream::app);
	//~ ReportFile << "height\tangle\tlambda\tmu\tabovemuBottom\tabove\n";
	bool piac = PointIsAboveCurve(argc,argv,h,a,study_num);
	ReportFile.close();

	//now increase height gradually; decrease angle if flow starts
	if (h<hMax) {
		if (piac) a-=.5;
		else h+=2;
		//write arg file
		std::stringstream command;
		command << "echo " << study_num << " " << h << " " << hMax << " " << a << " > arg";
		std::cout << command.str() << std::endl;
		return system (command.str().c_str());
	} else {
		return 0;	
	}
}

int main(int argc, char *argv[])
{
	if (argc<5) { 
	  std::cout << "Please specify Study, Height, max. Height and Angle" << std::endl; 
		exit(-1); 
	}
	HstopCurve(argc,argv,atoi(argv[1]),atof(argv[2]),atof(argv[3]),atof(argv[4]));	
}
