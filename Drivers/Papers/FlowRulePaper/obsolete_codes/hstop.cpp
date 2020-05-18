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
	problem.solve();
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

void HstopCurve(int study_num)
{
	Mdouble hStart = 4; //height at which bisection algorithm starts
	Mdouble aStart = 21; //angle at which bisection algorithm starts
	
	Mdouble hMax = 60; //height at which bisection algorithm stops
	Mdouble dh = 2; //how much I decrease the height after a point is found
	Mdouble da = .5; // how exact I want the angle to be determined
	Mdouble daMax = 1; //how much (in absolute) I increase the angle if the chute still stops

	stringstream name;
	name << "Report" << study_num;
	ReportFile.open(name.str().c_str(),fstream::out);
	ReportFile << "height\tangle\tlambda\tmu\tabovemuBottom\tabove\n";

	//we start at (aStart,hStart) and find the largest non-flowing angle for hStart
	//first we search in large steps daMax
	
	//find a flowing angle aMax and stopping angle aMin
	Mdouble aMax = aStart;
	while(!PointIsAboveCurve(hStart,aMax,study_num)) aMax+=daMax;
	Mdouble aMin=aMax-daMax;
	if (aMax==aStart) while(PointIsAboveCurve(hStart,aMin,study_num)) {aMin-=daMax; aMax-=daMax;}
	//bisect until aMin is highest stopping angle
	while (aMax-aMin>da*1.99) {
		Mdouble a=(aMin+aMax)/2;
		if (PointIsAboveCurve(hStart,a,study_num)) aMax=a;
		else aMin=a;
	}
	//now increase height gradually; decrease angle if flow starts
	Mdouble a = aMin;
	for (Mdouble h=hStart+dh; h<hMax; h+=dh) {
		while (PointIsAboveCurve(h,a,study_num)) a-=da;
	}
	
	ReportFile.close();
	return;	
}

int main(int argc, char *argv[])
{
	if (argc>1) HstopCurve(atoi(argv[1]));	
	else exit(-1);
}
