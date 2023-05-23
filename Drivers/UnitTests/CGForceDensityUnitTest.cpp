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

#include "Mercury3D.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "CG/CG.h"

/*! 
 * \brief Tests if the force density between two species has the correct sign assigned
 * \details The problem consists of two identical particles of different (but equal) species
 * in contact. 
 */
class CGForceDensityUnitTest : public Mercury3D
{
public:

    CGForceDensityUnitTest()
    {
        setName("CGForceDensityUnitTest");

        //set gravity and the species properties
        setGravity(Vec3D(-0.1, 0, -1.0));

        LinearViscoelasticSpecies species;
        species.setDensity(1.0);
        species.setCollisionTimeAndRestitutionCoefficient(.01, .1, 1.);

        speciesHandler.copyAndAddObject(species);
        speciesHandler.copyAndAddObject(species);
    }

    void setupInitialConditions() override {
        //define the domain

        setMax(Vec3D(1, 0.5, 0.5));
        setMin(Vec3D(-1, -0.5, -0.5));

        //define the two particles
        SphericalParticle p0;
        p0.setRadius(0.5);

        p0.setSpecies(speciesHandler.getObject(0));
        p0.setPosition(Vec3D(-0.49, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);
    
        p0.setSpecies(speciesHandler.getObject(1));
        p0.setPosition(Vec3D(0.49, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);
    }
};

int main(int argc, char *argv[])
{
    //declare the DPM problem and set the name
    CGForceDensityUnitTest problem;
    problem.setTimeMax(0.1); //a total of 5000 time steps
    problem.setTimeStep(0.1);
    problem.setFileType(FileType::NO_FILE);

    //define different coarse-graining objects (which should all result in the same mean values)
    auto cg0 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
    cg0->setTimeMin(0.1);
    cg0->selectSpecies(0);
    auto cg1 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
    cg1->setTimeMin(0.1);
    cg1->selectSpecies(1);
    
    //run the simulation
    problem.solve(argc, argv);

    Vec3D cg0ForceDensityResult = static_cast<CG<CGCoordinates::O,CGFunctions::Lucy,CGFields::StandardFields>*>(cg0)->evaluateAverage().getInteractionForceDensity();
    Vec3D cg1ForceDensityResult = static_cast<CG<CGCoordinates::O,CGFunctions::Lucy,CGFields::StandardFields>*>(cg1)->evaluateAverage().getInteractionForceDensity();

    logger(INFO, "ForceDensity species 0 %, forceDensity species 1 %, sum %", cg0ForceDensityResult, cg1ForceDensityResult, (cg0ForceDensityResult + cg1ForceDensityResult));
    logger.assert_always(cg0ForceDensityResult.getLength() > 10 && cg1ForceDensityResult.getLength() > 10, "evaluated force density is too small");
    logger.assert_always(cg0ForceDensityResult.X < -1 && cg1ForceDensityResult.X > 1, "force density does not point in the right direction");
    logger.assert_always( (cg0ForceDensityResult + cg1ForceDensityResult).getLength() < 1e-10, "evaluated force densities do not sum to 0");
}

