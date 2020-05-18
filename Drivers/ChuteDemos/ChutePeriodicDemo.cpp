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

#include <Boundaries/CubeInsertionBoundary.h>
#include <Boundaries/ChuteInsertionBoundary.h>
#include "Chute.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"

class ChutePeriodicDemo : public Chute
{
public:
    ChutePeriodicDemo()
    {
        setName("ChutePeriodicDemo");
        
        //Set and add the particle-species.
        LinearViscoelasticFrictionSpecies species;
        species.setDensity(6 / constants::pi);
        species.setCollisionTimeAndRestitutionCoefficient(5e-3, 0.88, 1);
        species.setSlidingDissipation(species.getDissipation());
        species.setSlidingStiffness(2. / 7. * species.getStiffness());
        species.setSlidingFrictionCoefficient(0.5);
        speciesHandler.copyAndAddObject(species);
        
        // Chute properties
        setChuteLength(30);
        setChuteWidth(10);
        setInflowHeight(10);
        setInflowParticleRadius(0.5);
        setFixedParticleRadius(getInflowParticleRadius());
        setRoughBottomType(MULTILAYER);
        
        //make chute periodic in y-direction:
        makeChutePeriodic();
        //make the chute periodic in x-direction:
        PeriodicBoundary b0;
        b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        
        //Set time-related properties
        setTimeMax(10);
        setTimeStep(1.0e-4); //(1/50th of the collision time)
        // Write 100 output files in total.
        setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(51, getTimeMax(), getTimeStep()));
    }
    
    ///add particles
    void setupInitialConditions() override
    {
        Chute::setupInitialConditions();
        addFlowParticlesCompactly();
        setParticlesWriteVTK(true);
    }
    
};

/** A rough inclined plane with periodic boundary conditions.
 **/
int main(int argc, char* argv[])
{
    //Print description
    logger(INFO, "\nDescription: A rough inclined plane with periodic boundary conditions.");
    
    // Problem parameters
    ChutePeriodicDemo problem;
    
    //solve
    problem.solve(argc, argv);
}
