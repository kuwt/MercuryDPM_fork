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

#include <iomanip>

#include "SilbertPeriodic.h"

class vibratedBed : public SilbertPeriodic {
public:
    
    vibratedBed() : SilbertPeriodic()
    {
        amplitude_=0.0;
        frequency_=0.0;
        setRoughBottomType(MULTILAYER);
    }
    
    void run(std::vector<Mdouble> studyNumber, int argc, char *argv[])
	{
		setInflowHeight(studyNumber[1]);
		setChuteAngle(studyNumber[2]);
        set_study(studyNumber[0]);
        std::stringstream ss;
        ss << "vibratedBed" << getName() << "A" << amplitude_ << "F" << frequency_;
        setName(ss.str());
        std::cout << "File name: " << getName() << std::endl;
        readArguments(argc, argv);
		writeRestartFile();
        solve();
	}
    
    void actionsAfterTimeStep() override {
        Vec3D velocity = getPrescribedVelocity();
        for (BaseParticle* p : particleHandler)
            if (p->isFixed())
                p->setVelocity(velocity);
        for (BaseWall* w : wallHandler)
            w->setVelocity(velocity);
    }

    Vec3D getPrescribedVelocity()
    {
        return Vec3D(0.0,0.0,amplitude_*2.0*constants::pi*frequency_*std::cos(2.0*constants::pi*frequency_*getTime()));
    }

    void setShakingStrength(Mdouble S)
    {
        //S=A*A*w*w
        frequency_ = std::sqrt(S/amplitude_/amplitude_);
    }

    void setFrequency(Mdouble frequency)
    {
        frequency_ = frequency;
    }

    void setAmplitude(Mdouble amplitude)
    {
        amplitude_ = amplitude;
    }

    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
            << ", tmax="  << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
            << ", N="     << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects()
            << ", theta=" << std::setprecision(3) << std::left << std::setw(6) << getChuteAngleDegrees()
            << ", Ene="   << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy() / getElasticEnergy()
            //~ << ", finish by " << setprecision(3) << left << setw(6) << timer.getFinishTime(t)
            << ". " << std::endl;
    }

private:
    Mdouble amplitude_;
    Mdouble frequency_;
};

int main(int argc, char *argv[])
{
	vibratedBed problem;
	//set case, height, angle to given or default values
	std::vector<Mdouble> studyNumber;
	studyNumber.resize(3);
    const unsigned int numArgs = 5;

	if (argc>numArgs)
    {
        studyNumber[0]=atof(argv[1]);
        studyNumber[1]=atof(argv[2]);
        studyNumber[2]=atof(argv[3]);
        problem.setAmplitude(atof(argv[4]));
        problem.setFrequency(atof(argv[5]));
        problem.setTimeStep(problem.getTimeStep()*2.0);
        problem.setSaveCount(50.0/problem.getTimeStep());
        problem.eneFile.setSaveCount(0.5/problem.getTimeStep());
        problem.setTimeMax(2000.0);
        problem.dataFile.setFileType(FileType::ONE_FILE);
        problem.restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
        problem.run(studyNumber,argc-numArgs,argv+numArgs);
    }
    else
    {
        std::cout << "Not enough input arguments given (./vibratedBed $study $height $angle $amplitude $frequency); " << std::endl
            << "using demo values (equivalent to ./vibratedBed 4 5 10 0.02 1 -tmax 1000 -dt 0.001)" << std::endl;
        studyNumber[0]=4;
        studyNumber[1]=0.1;
        studyNumber[2]=10;
        problem.setAmplitude(0.01);
        problem.setFrequency(1.0);
        problem.setRoughBottomType(MONOLAYER_DISORDERED);

        problem.setTimeStep(problem.getTimeStep()*10.0);
        problem.setTimeMax(2.0);
        problem.setSaveCount(1e1);
        problem.eneFile.setSaveCount(1e2);
        problem.dataFile.setFileType(FileType::ONE_FILE);
        problem.run(studyNumber,1,argv);
    }
}
