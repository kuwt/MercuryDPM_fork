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

#include "StatisticsVector.h"
#include "CG/TimeAveragedCG.h"
#include "BidispersedChute.h"

class ChuteForVelocityPaper : public BidispersedChute
{
public:
    ChuteForVelocityPaper(BidispersedChuteParameters parameters, Mdouble lengthX, Mdouble widthY) :
            BidispersedChute(parameters)
    {
        setChuteLength(lengthX);
        setChuteWidth(widthY);
        setMin(0, 0, 0);
        setMax(lengthX, widthY, 4);
        setRoughBottomType(RoughBottomType::MULTILAYER);
        std::stringstream name;
        name << 'H' << parameters.getInflowHeight()
             << "A" << parameters.getAngleInDegrees()
             << "Length" << lengthX
             << "Width" << widthY;
        setName(name.str());
        setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::NO_FILE);
        //dataFile.setSaveCount(10000);
        restartFile.setSaveCount(50./getTimeStep());
        eneFile.setSaveCount(0.1/getTimeStep());
        
        setTimeMax(2000);
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
    
    void actionsAfterSolve()
    {
        writeOutputFiles();
    }
};

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cout << "Please specify Height, size in x-direction, and size in y-direction" << std::endl;
        exit(-1);
    }
    Mdouble height = atof(argv[1]);
    auto parameters = BidispersedChuteParameters(height, 26, 0);
    parameters.setLargeParticleRadius(0.5);
    parameters.setSmallParticleRadius(0.5);
    parameters.setFixedParticleRadius(0.5);
    
    
    ChuteForVelocityPaper problem(parameters, atof(argv[2]), atof(argv[3]));
    problem.solve();
    return 0;
}
