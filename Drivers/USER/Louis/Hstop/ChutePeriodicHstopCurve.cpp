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

#include "ChutePeriodic.h"
//#include "ChutePeriodicHstopCurveParameters.h"

// todo merge the actionsBefore/AfterTimeStep into ChutePeriodic.h (e.g. bool checkFlowState, isFlowing, nextTimeToCheck),
// todo and specify time and file settings in runOnePeriodicChuteFlow() when a problem is created?
class ChutePeriodicHstopCurve : public ChutePeriodic
{
public:
    ChutePeriodicHstopCurve()
    {
        setName("ChutePeriodicHstopCurve");

        // necessary variables
        isFlowing = true;
        nextTimeToCheck = 10;

        // time related settings
        setTimeMax(500);
        setTimeStep(1.0e-4); //(1/50th of the collision time)
        setSaveCount(getTimeMax()/getTimeStep());
        eneFile.setSaveCount(0.1/getTimeStep());

        // file related settings
        fStatFile.setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        eneFile.setFileType(FileType::ONE_FILE);

        // other settings
        setRoughBottomType(RoughBottomType::MONOLAYER_TRIANGULAR);
    }

    void actionsBeforeTimeStep() override
    {
        // inherit from "parent function"
        ChutePeriodic::actionsBeforeTimeStep();

        // report current time info to screen at every time unit
        if (ceil(getTime()) != ceil(getTime() + getTimeStep()))
            printTime();
    }

    void actionsAfterTimeStep() override
    {
        // inherit from "parent function"
        ChutePeriodic::actionsAfterTimeStep();

        // check flow state every 10 time units
        if (getTime() > nextTimeToCheck)
        {
            // update flow state and print out
            if (getKineticEnergy() / getElasticEnergy() < 1e-5)
                isFlowing = false;
            std::cout << "ene_kin/ene_ele=" << getKineticEnergy() / getElasticEnergy()
                      << ", isFlowing = " << isFlowing << std::endl;

            // update nextTimeToCheck
            nextTimeToCheck += 10;
        }
    }

    bool continueSolve() const override
    {
        return (isFlowing);
    }

    void solve_analytic()
    {
        setRoughBottomType(FLAT);
        setTimeMax(1e-4);
        solve();
        isFlowing = (tan(getChuteAngle()) < .2 + .2 * exp(getZMax()));
        isFlowing = (getZMax() > 80 - 4 * getChuteAngleDegrees());
        isFlowing = (getInflowHeight() > 80 - 4 * getChuteAngleDegrees());
    };

    bool getIsFlowing() {
        return isFlowing;
    }

private:
    bool isFlowing;
    Mdouble nextTimeToCheck;
};

// ------------------------------------

class HstopCurve
{
public:

    HstopCurve()
    {
        // three main parameters
        angleFirst_ = 26;
        heightFirst_ = 40;
        heightSecond_ = 4;

        // other variables
        isTest_ = false;
        pointIsAboveCurve_ = true;
        a = angleFirst_;
        h = std::min(heightFirst_,heightSecond_);
        hMax = std::max(heightFirst_,heightSecond_);

        // other parameters
        fixParticleRadius_ = 0.5;
    }

    HstopCurve(Mdouble angleFirst, Mdouble heightFirst, Mdouble heightSecond)
    {
        // three main parameters
        angleFirst_ = angleFirst;
        heightFirst_ = heightFirst;
        heightSecond_ = heightSecond;

        // other variables
        // todo: repeated with the first constructor?
        isTest_ = false;
        pointIsAboveCurve_ = true;
        a = angleFirst_;
        h = std::min(heightFirst_,heightSecond_);
        hMax = std::max(heightFirst_,heightSecond_);

        // other parameters
        fixParticleRadius_ = 0.5;
    }

    void search()
    {
        //initialize search routine
        initializeSearchRoutine();

        //search loop
        while (h < hMax)
        {
            //run one case and check flow state
            runOnePeriodicChuteFlow();

            //write result into a report
            writeToHstopCurveReport();

            //update h and a
            updateSearchRoutine();
        }
    }

    // setters (options for hstop curve searching including bottom properties etc.)

    void setIsTest(bool isTest)
    {
        isTest_ = isTest;
    }

//    void passArguments(int argc, char *argv[])
//    {
//        argc_ = argc - 3;
//        argv_ = argv + 3;
//    }

private:

    void initializeSearchRoutine()
    {
        // initialize increments
        // initialize a, h, hMax
    }

    void updateSearchRoutine()
    {
        if (pointIsAboveCurve_)
        {
            a -= 0.5;
        }
        else
        {
            if (h < 6 || a > 24)
                h *= 1.1;
            else
                h += 2;
        }
        if (h >= 6)
            h = round(h);
    }

    void runOnePeriodicChuteFlow()
    {
        //name the case and print out
        std::stringstream name;
        name << "H" << h << "A" << a;
        std::stringstream logInfo;
        logInfo << "\n----------------------------------\n"
                << "running one periodic chute flow: " << name.str().c_str() << std::endl;
        logger(INFO, logInfo.str().c_str());

        //create the case with h and a
        ChutePeriodicHstopCurve problem;
        problem.setName(name.str().c_str());
        problem.setInflowHeight(h);
        problem.setChuteAngle(a);

        double h_min = 10; //10 particle diameters
        problem.setZMax(std::min(h_min,1.2*h));
        problem.setNumberOfDomains({1,2,2}); //12 cores

        // todo parameter study is specified here
        // todo should pass argc and argv to here
        //parameters other than h and a, such as bottom properties
        problem.setFixedParticleRadius(fixParticleRadius_);
//        problem.readArguments(argc_,argv_);
        // run
        if (isTest_)
        {
            problem.solve_analytic();
        }
        else
        {
            problem.solve();
        }

        //update pointIsAboveCurve_
        pointIsAboveCurve_ = problem.getIsFlowing();
    }

    void writeToHstopCurveReport()
    {
        //name the report file
        std::stringstream filename;
        filename << "Report"; // followed by major parameters, e.g. L1 or L1M05B05 or Ra05

        //write (append) to the file
        std::fstream ReportFile;
        ReportFile.open(filename.str().c_str(), std::fstream::out | std::fstream::app);
        ReportFile << h << "\t" << a << "\t" << pointIsAboveCurve_ << std::endl;
        ReportFile.close();
    }

    // variables that do not change during run
    bool isTest_;
    Mdouble angleFirst_;
    Mdouble heightFirst_;
    Mdouble heightSecond_;

    // variables that control searching loop
    bool pointIsAboveCurve_;
    Mdouble hMax, h, a;

    // variables to be passed to each problem (except hMax, h, a)
    // todo alternative: define varibles that receive argc and argv?
    Mdouble fixParticleRadius_;
//    int argc_;
//    char *argv_[];

};

//====================================================

int main(int argc, char* argv[])
{
    // todo how to pass options from arguments?
    HstopCurve hstopCurve(20,4,40);
//    hstopCurve.passArguments(argc,argv);
    hstopCurve.setIsTest(false);
    hstopCurve.search();
    return 0;
}