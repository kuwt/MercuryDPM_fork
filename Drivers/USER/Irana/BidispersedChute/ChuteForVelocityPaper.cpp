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

#include "StatisticsVector.h"
#include "CG/TimeAveragedCG.h"
#include "BidispersedChute.h"

class ChuteForVelocityPaper : public BidispersedChute
{
public:
    ChuteForVelocityPaper(BidispersedChuteParameters parameters) : BidispersedChute(parameters)
    {
        setChuteLength(20);
        setChuteWidth(10);
        setMin(0, 0, 0);
        setMax(20, 10, 4);
        setFileType(FileType::ONE_FILE);
        setRoughBottomType(RoughBottomType::MULTILAYER);
        std::stringstream name;
        name << 'H' << parameters.getInflowHeight()
             << "A" << parameters.getAngleInDegrees()
             << "Sr" << (int) std::pow(parameters.getLargeParticleRadius() / parameters.getSmallParticleRadius(), 3);
        setName(name.str());
        setSaveCount(50./getTimeStep());
        if (parameters.getConcentrationSmall() < 1e-10)
        {
            setTimeMax(500);
        }
        else
        {
            setTimeMax(2000);
        }
        std::string filename = "COM"
                               + std::to_string(parameters.getInflowHeight()) + "_"
                               + std::to_string((int) std::pow(parameters.getLargeParticleRadius() /
                                                               parameters.getSmallParticleRadius(), 3) )+ "_"
                               + std::to_string(parameters.getAngleInDegrees());
        comFile.open(filename, std::ofstream::out);
        comFile << std::setw(12);
        comFile << "time" << std::setw(12) << "com_flow" << std::setw(12) << "com_large" << std::setw(12)
                <<"com_small" << std::endl;
    }
    
    void setupInitialConditions() override
    {
        BidispersedChute::setupInitialConditions();
        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed())
            {
                p->setSpecies(speciesHandler.getObject(0));
            }
        }
        logger(INFO, "Name: %, savecount restart file: %", getName(), restartFile.getSaveCount());
    }
    
    ~ChuteForVelocityPaper()
    {
        comFile.close();
    }
    
    void actionsAfterSolve()
    {
        writeOutputFiles();
    }
    
    ///compute the centre of mass of a) all flow particles b) all large particles c) all small particles
    ///also check if the flow is arrested.
    ///note, that we do not check every time step, but only every 1t.
    void actionsAfterTimeStep() override
    {
        static Mdouble nextCheckedTime = 5;
        if (getTime() > nextCheckedTime)
        {
            //compute COM of species 1
            Mdouble com1 = 0;
            int n1 = 0;
            //compute COM of species 2
            Mdouble com2 = 0;
            int n2 = 0;
            //compute COM of species 1 and 2
            Mdouble com = 0;
            int n = 0;
            for (const BaseParticle* const p : particleHandler)
            {
                if (p->getSpecies()->getIndex() == 1)
                {
                    com1 += p->getPosition().Z;
                    n1++;
                    com += p->getPosition().Z;
                    n++;
                }
                else if (p->getSpecies()->getIndex() == 2)
                {
                    com2 += p->getPosition().Z;
                    n2++;
                    com += p->getPosition().Z;
                    n++;
                }
            }
            com /= n;
            com1 /= n1;
            if (n2 > 0)
                com2 /= n2;
            nextCheckedTime += 5;
            logger(INFO, "COM species 1: %, COM species 2: %", com1, com2);
            int width = 10;
            comFile << std::setw(12);
            comFile << getTime() << std::setw(12) << com << std::setw(12) << com1 << std::setw(12) << com2 << std::endl;
            if (getKineticEnergy() < 1e-5)
                logger(ERROR, "The flow has arrested");
        }
    }

private:
    
    std::ofstream comFile;
};

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cout << "Please specify Height, Angle, and mass-ratio" << std::endl;
        exit(-1);
    }
    Mdouble height = atof(argv[1]);
    Mdouble angle = atof(argv[2]);
    Mdouble sizeRatio = atof(argv[3]);
    auto parameters = BidispersedChuteParameters(height, angle, 0.5);
    parameters.setFixedParticleRadius(0.5);
    if (mathsFunc::isEqual(sizeRatio, 1, 1e-2))
    {
        parameters = BidispersedChuteParameters(height, angle, 0);
        parameters.setLargeParticleRadius(0.5);
        parameters.setSmallParticleRadius(0.5);
    }
    else
    {
        parameters.setSmallParticleRadius(1./(1 + std::cbrt(sizeRatio)));
        parameters.setLargeParticleRadius(1 - parameters.getSmallParticleRadius());
    }
    
    
    logger(INFO, "radii of particles: %, %", parameters.getSmallParticleRadius(), parameters.getLargeParticleRadius());
    ChuteForVelocityPaper problem(parameters);
    problem.solve();
    return 0;
}
