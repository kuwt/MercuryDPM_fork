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

#include "../BidispersedChute/BidispersedChute.h"

class SilbertHstop : public BidispersedChute
{
public:
    
    SilbertHstop(const BidispersedChuteParameters& bidispersedChuteParameters = BidispersedChuteParameters()) :
            BidispersedChute(bidispersedChuteParameters)
    {
        flowIsStillMoving = true;
        //these simulations take up to a day, so no restart files needed (otherwise generates a lot of data)
        setTimeMax(20);
        setSaveCount((unsigned int) (10*getTimeMax() / getTimeStep()));
        fStatFile.setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::NO_FILE);
        eneFile.setFileType(FileType::NO_FILE);
        statFile.setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::MULTIPLE_FILES);
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
    
    void resetAfterFlowing(Mdouble angleDecrement)
    {
        setChuteAngle(getChuteAngleDegrees() - angleDecrement);
        parameters.setAngle(getChuteAngle());
        setTime(0);
        setRestarted(true);
        logger(INFO, "New angle: %", getChuteAngleDegrees());
    }
    
    void resetAfterArresting(Mdouble extraHeight)
    {
        unsigned int extraLarge = parameters.getNumberOfLargeParticlesExtraAndUpdate(getChuteLength() * getChuteWidth(), extraHeight);
        unsigned int extraSmall = parameters.getNumberOfSmallParticlesExtraAndUpdate(getChuteLength() * getChuteWidth(), extraHeight);
        insertParticles(extraLarge, extraSmall);
        setInflowHeight(parameters.getInflowHeight());
        for (BaseParticle* p : particleHandler)
        {
            p->setVelocity(p->getVelocity() + Vec3D(.1,0,0));
        }
        setTime(0);
        flowIsStillMoving = true;
        setRestarted(true);
        logger(INFO, "new inflow height: %", getInflowHeight());
        logger(INFO, "Total number of particles: %", particleHandler.getNumberOfObjects());
    }
    
    ///\brief indicator of whether the flow is moving, or arrested.
    ///\details indicator of whether the flow is still moving, or arrested. Is always true for t < 10, afterwards it is
    /// true if the flow is flowing, false if the flow is arrested (kineticEnergy/elasticEnergy < 1e-5).
    bool flowIsStillMoving;
};

bool PointIsAboveCurve(SilbertHstop& problem, bool test)
{
    if (test)
    {
        //for test runs:
        logger(INFO, "testing\n\n\n\n\n");
        problem.testWithAnalyticFunction(problem.getInflowHeight(), problem.getChuteAngleDegrees());
    }
    else
    {
        //for real runs:
        problem.solve();
    }
    return false;
    return problem.flowIsStillMoving;
}

void hstopCurve(Mdouble h, Mdouble hMax, Mdouble angle, Mdouble phi, Mdouble sizeRatio)
{
    bool test = false;
    const Mdouble angleDecrement = 0.5;
    const Mdouble heightIncrement = 2;
    
    std::stringstream fileName;
    fileName << "Report_" << phi << "_" << sizeRatio;
    std::ofstream reportFile;
    reportFile.open(fileName.str().c_str(), std::fstream::out | std::fstream::app);
    BidispersedChuteParameters bidispersedChuteParameters(h, angle, phi);
    SilbertHstop problem(bidispersedChuteParameters);
    while (h < hMax)
    {
        logger(INFO, "check if point is above curve for h = %, a = % (phi = %, s = %)", h, angle, phi);
        std::stringstream name;
        name << "height_" << h
             << "_angle_" << angle
             << "_phi_" << phi
             << "_sizeRatio_" << sizeRatio;
        problem.setName(name.str().c_str());
        logger(INFO, "starting %", name.str());
    
        //now increase height gradually; decrease angle if flow starts
        bool flowing = PointIsAboveCurve(problem, test);
        if (flowing)
        {
            logger(INFO, "Flow keeps flowing");
            angle -= angleDecrement;
            problem.resetAfterFlowing(angleDecrement);
        }
        else
        {
            logger(INFO, "Flow arrests");
            h += heightIncrement;
            problem.resetAfterArresting(heightIncrement);
        }
        
        reportFile << name.str()
                   << "\t" << flowing
                   << std::endl;
    }
    reportFile.close();
}

int main(int argc, char* argv[])
{
    const Mdouble h = 4;
    const Mdouble hMax = 5;
    const Mdouble startingAngle = 24;
    const Mdouble concentrationSmall = 0.0;
    const Mdouble sizeRatio = std::cbrt(2);
    hstopCurve(h, hMax, startingAngle, concentrationSmall, sizeRatio);
    return 0;
}
