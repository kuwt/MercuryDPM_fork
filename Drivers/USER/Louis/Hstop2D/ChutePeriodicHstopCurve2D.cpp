//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include "ChutePeriodic2D.h"

class HstopCurve
{
public:

    HstopCurve()
    {
        a = 15;
        h = 5;
        hMax = 100;
        bt = ORDERED;
        srb = 1;
        eps = 0;
        pointIsAboveCurve_ = true;
    }

    HstopCurve(Mdouble angle, Mdouble heightMin, Mdouble heightMax)
    {
        a = angle;
        h = heightMin;
        hMax = heightMax;
        bt = ORDERED;
        srb = 1;
        eps = 0;
        pointIsAboveCurve_ = true;
    }

    void search()
    {

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

    // getters



    // setters

    void setBaseType(BottomType baseType)
    {
        bt = baseType;
    }

    void setFixedParticleRadiusRatio(Mdouble fixedParticleRadiusRatio)
    {
        srb = fixedParticleRadiusRatio;
    }

    void setFixedParticleSpacing(Mdouble fixedParticleSpacing)
    {
        eps = fixedParticleSpacing;
    }

private:

    void updateSearchRoutine()
    {
        if (pointIsAboveCurve_)
        {
            a *= 0.98;
        }
        else
        {
            h *= 1.2;
        }
    }

    void runOnePeriodicChuteFlow()
    {
        //name the case and print out
        std::stringstream name;
        name << std::fixed;
        name << "H" << std::setprecision(2) << h << "A" << std::setprecision(2) << a;
        std::stringstream logInfo;
        logInfo << "\n----------------------------------\n"
                << "running one periodic chute flow: " << name.str().c_str() << std::endl;
        logger(INFO, logInfo.str().c_str());

        //create the case with h and a
        ChutePeriodic2D problem;
        problem.setName(name.str().c_str());
        problem.setInflowHeight(h);
        problem.setChuteAngle(a);
        problem.setBottomType(bt);
        problem.setFixedParticleRadiusRatio(srb);
        problem.setFixedParticleSpacing(eps);
        problem.setTimeMax(200);
        problem.cgHandler.removeObject(1);
        problem.cgHandler.removeObject(0);
        problem.eneFile.setFileType(FileType::ONE_FILE);
        problem.dataFile.setFileType(FileType::NO_FILE);
        problem.fStatFile.setFileType(FileType::NO_FILE);
        problem.restartFile.setFileType(FileType::NO_FILE);

        double h_min = 10; //10 particle diameters
        problem.setZMax(std::min(h_min,1.2*h));

        problem.solve();

        //update pointIsAboveCurve_
        pointIsAboveCurve_ = problem.getIsFlowing();
    }

    void writeToHstopCurveReport()
    {
        //name the report file
        std::string btShort;
        switch(bt) {
            case FLAT : btShort = "F"; break;
            case ORDERED : btShort = "O"; break;
            case DISORDERED : default:  btShort = "D"; break;
        }
        std::stringstream filename;
        filename << "Hstop_" << btShort
                 << "_R" << srb
                 << "E" << eps;

        //write (append) to the file
        std::fstream ReportFile;
        ReportFile.open(filename.str().c_str(), std::fstream::out | std::fstream::app);
        ReportFile << h << "\t" << a << "\t" << pointIsAboveCurve_ << std::endl;
        ReportFile.close();
    }

    // variables that control searching loop
    bool pointIsAboveCurve_;
    Mdouble hMax, h, a;
    // variables for bottom
    BottomType bt;
    Mdouble srb, eps;

    // variables to be passed to each problem (except hMax, h, a)

};

//====================================================

int main(int argc, char* argv[])
{
    // read arguments
    if (argc < 6)
    {
        std::cout << "Please specify Height, HeightMax, Angle, BottomType, SizeRatio, and Spacing." << std::endl;
        exit(-1);
    }

    Mdouble h=strtod(argv[1],nullptr);
    Mdouble hmax=strtod(argv[2],nullptr);
    Mdouble a=strtod(argv[3],nullptr);
    auto bt= static_cast<int>(strtod(argv[4], nullptr));
    Mdouble srb=strtod(argv[5],nullptr);
    Mdouble eps=strtod(argv[6],nullptr);
    std::string btShort;
    BottomType bottomType;
    switch(bt) {
        case 0 : bottomType = FLAT; break;
        case 1 : bottomType = ORDERED; break;
        case 2 : default: bottomType = DISORDERED; break;
    }

    HstopCurve hstopCurve(a,h,hmax);
    hstopCurve.setBaseType(bottomType);
    hstopCurve.setFixedParticleRadiusRatio(srb);
    hstopCurve.setFixedParticleSpacing(eps);
    hstopCurve.search();

    return 0;
}