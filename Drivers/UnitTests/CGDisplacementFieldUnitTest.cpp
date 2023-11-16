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
#include "CG/TimeAveragedCG.h"
#include "CG/Fields/DisplacementField.h"

/*! 
 * \brief Tests if the displacement momentum density is correctly calculated.
 * \details The problem consists of one particle moving at constant velocity.
 */
class CGDisplacementFieldUnitTest : public Mercury3D
{
public:

    CGDisplacementFieldUnitTest()
    {
        setName("CGDisplacementFieldUnitTest");

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

        //define a moving particle
        SphericalParticle p0;

        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.1);
        p0.setPosition(Vec3D(-0.9, 0.0, 0.0));
        p0.setVelocity(Vec3D(0.9, 0.0, 0.0));

        BaseParticle* p1 = particleHandler.copyAndAddObject(p0);
    }
};

int main(int argc, char *argv[])
{
    //declare the DPM problem and set the name
    CGDisplacementFieldUnitTest problem;
    problem.setTimeMax(1);
    problem.setTimeStep(0.1);
    problem.setFileType(FileType::NO_FILE);

    //define different coarse-graining objects (which should all result in the same mean values)
    auto cg0 = problem.cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::O,CGFunctions::Lucy,CGFields::DisplacementField>());

    // Initially, the previous position is not set, so the displacement would not make sense at time step 0
    cg0->setTimeMin(0.2);
    cg0->selectSpecies(0);
    
    //run the simulation
    problem.solve(argc, argv);

    Vec3D displacementMomentum = cg0->evaluateAverage().getDisplacementMomentum();

    // Calculate the expected displacement momentum
    BaseParticle* p0 = problem.particleHandler.getObject(0);
    Vec3D diff = problem.getMax() - problem.getMin();
    Mdouble simulationVolume = diff.X * diff.Y * diff.Z;
    Vec3D expectedDisplacementMomentum = p0->getMass() * p0->getVelocity() / simulationVolume;

    logger.assert_always((displacementMomentum-expectedDisplacementMomentum).getLength() < 1e-10, "Displacement momentum is not correct");
}

