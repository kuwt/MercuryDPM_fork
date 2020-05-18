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

#include "GlasPeriodic.h"
using namespace std;

class SilbertHstop : public SilbertPeriodic 
{
public:

	SilbertHstop() {
		pointIsAboveCurve = true;
		//~ setTimeStep(getTimeStep()/3.);
		setSaveCount(1/getTimeStep()); 
		setTimeMax(500);
		fStatFile.setFileType(FileType::NO_FILE);
		dataFile.setFileType(FileType::NO_FILE);
	}

	Mdouble pointIsAboveCurve;
	
	void writeToEne()
	{
		Mdouble ene_kin = 0, ene_rot = 0, ene_gra = 0, mass_sum= 0, x_masslength=0, y_masslength=0, z_masslength=0;
		for (unsigned int i=0;i<particleHandler.getNumberOfObjects();i++) if (!particleHandler.getObject(i)->isFixed())
		{
			
			ene_kin += .5 * particleHandler.getObject(i)->getMass() * particleHandler.getObject(i)->getVelocity().getLengthSquared();
			ene_rot += particleHandler.getObject(i)->getRotationalEnergy();
			ene_gra -= particleHandler.getObject(i)->getMass() * Vec3D::dot(getGravity(),particleHandler.getObject(i)->getPosition());
			mass_sum +=particleHandler.getObject(i)->getMass();
			x_masslength +=particleHandler.getObject(i)->getMass()*particleHandler.getObject(i)->getPosition().X;
			y_masslength +=particleHandler.getObject(i)->getMass()*particleHandler.getObject(i)->getPosition().Y;
			z_masslength +=particleHandler.getObject(i)->getMass()*particleHandler.getObject(i)->getPosition().Z;
		} //end for loop over Particles

		///todo{Why is there a +6 here?}
		static int width = eneFile.getFstream().precision() + 6;
		eneFile.getFstream()  << setw(width) << getTime()
			<< " " << setw(width) << ene_gra
			<< " " << setw(width) << ene_kin
			<< " " << setw(width) << ene_rot
			<< " " << setw(width) << getElasticEnergy()
			<< " " << setw(width) << (mass_sum?x_masslength/mass_sum:NAN)
			<< " " << setw(width) << (mass_sum?y_masslength/mass_sum:NAN) 
			<< " " << setw(width) << (mass_sum?z_masslength/mass_sum:NAN)
			<< endl;
		
		cout << "ene_kin/ene_ele=" << ene_kin/getElasticEnergy() << ", piac=" << pointIsAboveCurve << endl;
		if (ene_kin/getElasticEnergy()<1e-5 && getTime()>10) pointIsAboveCurve=false;
		
		resetElasticEnergy();
	}

	bool continueSolve() {
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

fstream ReportFile;

bool PointIsAboveCurve(Mdouble h, Mdouble a, int study_num) {
	cout << "PointIsAboveCurve(" << study_num << ", h=" << h << ", a=" << a << ")" << endl;
	SilbertHstop problem;
	problem.setInflowHeight(h) ;
	problem.setChuteAngle(a);
	problem.set_study(study_num);
	stringstream name;
	name << "H" << problem.getInflowHeight() 
		 << "A" << problem.getChuteAngleDegrees() 
		 << "L" << round(100.*problem.getFixedParticleRadius()*2.)/100.
		 << "M" << problem.getSlidingFrictionCoefficient()
		 << "B" << problem.getSlidingFrictionCoefficientBottom();
	problem.setName(name.str().c_str());
	cout << "starting " << name.str().c_str() << endl;
	// if you want more restart output, use these options
	//problem.restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
	// for real runs:
	//~ problem.set_Hertzian(true);
	problem.solve();
	// for test runs:
	//problem.solve_analytic();
	stringstream com;
	ReportFile   << problem.getInflowHeight() 
		 << "\t" << problem.getChuteAngleDegrees() 
		 << "\t" << round(100.*problem.getFixedParticleRadius()*2.)/100.
		 << "\t" << problem.getSlidingFrictionCoefficient()
		 << "\t" << problem.getSlidingFrictionCoefficientBottom()
		 << "\t" << problem.pointIsAboveCurve
		 << endl;
	return problem.pointIsAboveCurve;
}

int HstopCurve(string cmd UNUSED, int study_num, Mdouble h, Mdouble hMax, Mdouble a)
{
	cout << "restart at study_num=" << study_num << ", h=" << h << ", a=" << a << endl;
	
	stringstream name;
	name << "Report" << study_num;
	ReportFile.open(name.str().c_str(),fstream::out | fstream::app);
	//~ ReportFile << "height\tangle\tlambda\tmu\tabovemuBottom\tabove\n";
	bool piac = PointIsAboveCurve(h,a,study_num);
	ReportFile.close();

	//now increase height gradually; decrease angle if flow starts
	if (h<hMax) {
		if (piac) a-=.5;
		else h+=2;
		stringstream command;
		// for msm1:
		command << "nohup nice -18 ../hstopGlass_StudyHeightHmaxAngle.exe " << study_num << " " << h << " " << hMax << " " << a << " &>>outS" << study_num << "H" << h << "HMax" << hMax << "A" << a << " &";
		// for einder:
		//command << "name="<<cmd<<" && cd ${name/hstop_StudyHeightHmaxAngle.exe/} && ~/clusterscriptexecute $name " << study_num << " " << h << " " << hMax << " " << a << "&";
		cout << command.str() << endl;
		return system (command.str().c_str());
	} else {
		return 0;	
	}
}

int main(int argc, char *argv[])
{
	if (argc<4) { cout << "Please specify Study, Height, max. Height and Angle" << endl; exit(-1); }
	HstopCurve(argv[0],atoi(argv[1]),atof(argv[2]),atof(argv[3]),atof(argv[4]));	
}
