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

/** The file hstop_StudyHeightHmaxAngle_new.cpp defines a simple bisection algorithm to determine
 * the hstop curve, i.e. at which angle and height a steady chute flow is not possible anymore and the flow stops.
 * It starts at a low height and determines the highest angle at which the chute flow subsides,
 * then increases height and repeats the process (it takes into account that the hstop(theta) curve is monotonically
 * decreasing). It does this by calling PointIsAboveCurve, which runs a simulation of a fixed height and angle
 * and returns 0 is the flow subsides before tmax and 1 if not
 **/

#include <algorithm>
#include "../BidispersedChute/BidispersedChute.h"

class HStopFixedHeight : public BidispersedChute
{
public:
    
    HStopFixedHeight(const BidispersedChuteParameters& bidispersedChuteParameters = BidispersedChuteParameters()) :
            BidispersedChute(bidispersedChuteParameters)
    {
        flowIsStillMoving = true;
        //these simulations take up to a day, so no restart files needed (otherwise generates a lot of data)
        setTimeMax(1000);
        fStatFile.setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::NO_FILE);
        eneFile.setFileType(FileType::ONE_FILE);
        statFile.setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::NO_FILE);
    }
    
    
    void actionsAfterTimeStep() override
    {
        BidispersedChute::actionsAfterTimeStep();
        if (getKineticEnergy()  < 1e-5 * getElasticEnergy() && getTime() > 10)
            flowIsStillMoving = false;
    }
    
    bool continueSolve() const override
    {
        return flowIsStillMoving;
    }
    
    void testWithAnalyticFunction(Mdouble h, Mdouble a)
    {
        Mdouble func = 20 * (tan(a / 180. * constants::pi) - tan(30. / 180. * constants::pi)) /
                       (tan(20. / 180. * constants::pi) - tan(30. / 180. * constants::pi));
        flowIsStillMoving = (h > func);
    }
    
    void changeAngle(Mdouble angleNew)
    {
        setChuteAngle(angleNew);
        parameters.setAngle(getChuteAngle());
        logger(INFO, "New angle: %", getChuteAngleDegrees());
    }
    
    void actionsAfterSolve() override
    {
        writeRestartFile();
    }
    
    ///\brief indicator of whether the flow is moving, or arrested.
    ///\details indicator of whether the flow is still moving, or arrested. Is always true for t < 10, afterwards it is
    /// true if the flow is flowing, false if the flow is arrested (kineticEnergy/elasticEnergy < 1e-5).
    bool flowIsStillMoving;
};

void findCriticalAngle(Mdouble h, Mdouble concentrationSmall, Mdouble sizeRatio, Mdouble largeParticleRadius)
{
    const bool test = false;
    Mdouble angleDecrement = 1;
    const Mdouble minAngleDecrement = .1;
    Mdouble currentAngle = 25;
    Mdouble minFlowingAngle = currentAngle;
    std::string minFlowingName = "";
    std::vector<Mdouble> arrestedAngles;
    while (angleDecrement > minAngleDecrement)
    {
        BidispersedChuteParameters parameters = BidispersedChuteParameters(h, currentAngle, concentrationSmall);
        parameters.setLargeParticleRadius(largeParticleRadius);
        HStopFixedHeight problem(parameters);
        std::stringstream name;
        logger(INFO, "check if point is above curve for h = %, a = % (phi = %, s = %)", h, currentAngle, concentrationSmall, sizeRatio);
        name << "height_" << h
             << "_angle_" << currentAngle
             << "_phi_" << concentrationSmall
             << "_sizeRatio_" << sizeRatio
             << "_largeParticleRadius_" << largeParticleRadius;
        logger(INFO, "starting %", name.str());
    
        if (minFlowingName.size() > 0)
        {
            problem.readRestartFile(minFlowingName);
            problem.setTime(0);
        }
        problem.setName(name.str());
        if (test)
        {
            problem.testWithAnalyticFunction(h, currentAngle);
        }
        else
        {
            problem.changeAngle(currentAngle);
            problem.solve();
        }
        if (problem.flowIsStillMoving)
        {
            minFlowingAngle = currentAngle;
            currentAngle -= angleDecrement;
            minFlowingName = name.str();
            if (std::find(arrestedAngles.begin(), arrestedAngles.end(),currentAngle) != arrestedAngles.end())
            {
                angleDecrement /= 2;
                currentAngle = minFlowingAngle - angleDecrement;
            }
            logger(INFO, "flow flows, new angle %\n\n", currentAngle);
        }
        else
        {
            arrestedAngles.push_back(currentAngle);
            angleDecrement /= 2;
            currentAngle = minFlowingAngle - angleDecrement;
            logger(INFO, "flow arrests, new angle %\n\n", currentAngle);
        }
    }
    std::stringstream fileName;
    fileName << "Report_" << concentrationSmall << "_" << sizeRatio;
    std::ofstream reportFile;
    reportFile.open(fileName.str().c_str(), std::fstream::out | std::fstream::app);
    reportFile << "height \t" << h
               << "\tangle \t" << minFlowingAngle
               << "\tphi \t" << concentrationSmall
               << "\tsizeRatio\t" << sizeRatio
               << "\tlargeParticleRadius\t" << largeParticleRadius << "\n";
}

int main(int argc, char* argv[])
{
    logger.assert_always(argc > 2, "Please specify height, large particle radius");
    const Mdouble h = atof(argv[1]);
    const Mdouble largeParticleRadius = atof(argv[2]);
    const Mdouble concentrationSmall = 0.0;
    const Mdouble sizeRatio = std::cbrt(1);
    findCriticalAngle(h, concentrationSmall, sizeRatio, largeParticleRadius);
    return 0;
}
