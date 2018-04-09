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

#include "Mercury3D.h"
#include "Math/ExtendedMath.h"
#include "Species/HertzianViscoelasticMindlinSpecies.h"
#include "CG/CG.h"

class CreateDataAndFStatFiles : public DPMBase
{
public:
    void setupInitialConditions() override
    {
        setName("CGReadFiles");
        dataFile.setFileType(FileType::MULTIPLE_FILES);
        fStatFile.setFileType(FileType::MULTIPLE_FILES);
        setSaveCount(10);

        HertzianViscoelasticMindlinSpecies species;
        species.setDensity(6.0 / constants::pi);
        speciesHandler.copyAndAddObject(species);

        BaseParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(0.5);
        p.setPosition(Vec3D(0.0, 0.0, 0.0));
        p.setVelocity(Vec3D(-1.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p);

        setTimeStep(0.01);
        setTimeMax(100*getTimeStep());
        setMin(-0.5, -0.5, -0.5);
        setMax(0.5, 0.5, 0.5);
    }
};


class CGReadFiles : public DPMBase
{
public:

    void test()
    {
        setName("CGReadFiles2");
        auto cg = cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
        cgHandler.restartAndEvaluateDataFiles("CGReadFiles");
        Mdouble density = cg->getPoint(0).getDensity();
        std::cout << getTime();
        logger.assert_always(mathsFunc::isEqual(density,1,1e-10),
                "Wrong density: % (should be 1)",density);
        logger(INFO,"Test successful");
    }
};

int main()
{
    logger(INFO,"Testing if data and fstat files are successfully written, and read by the cgHandler.");

    CreateDataAndFStatFiles problem0;
    problem0.solve();

    CGReadFiles problem1;
    problem1.test();
    return 0;
}
