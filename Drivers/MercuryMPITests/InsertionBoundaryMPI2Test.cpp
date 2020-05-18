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

#include <Particles/LiquidFilmParticle.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Species/LinearViscoelasticSpecies.h>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"
using constants::pi;
using mathsFunc::cubic;

class InsertionBoundaryMPI2Test : public Mercury3D {
public:

    //sets up name, domain, parallel decomposition, time step, repulsive-adhesive species
    InsertionBoundaryMPI2Test() {
        setName("InsertionBoundaryMPI2Test");
        setMin(-Vec3D(5,1,1));
        setMax(Vec3D(5,1,1));
        setNumberOfDomains({NUMBER_OF_PROCESSORS,1,1});
        setTimeStep(1e-4);
        setTimeMax(1.0);

        LinearViscoelasticSpecies species;
        species.setDensity(6.0/pi);
//        species.setPlasticParameters(2e5,10e5,2e5,.5);
        species.setStiffness(2e5);
        species.setDissipation(25);
//        species.setSlidingFrictionCoefficient(0.5);
//        species.setSlidingDissipation(2./7.*species.getDissipation());
//        species.setSlidingStiffness(2./7.*species.getStiffness());
//        species.setLiquidBridgeVolumeMax(1e-10);
//        species.setContactAngle(0);
//        double effectiveRadius = 0.25;
//        species.setSurfaceTension(1./(2.0 * constants::pi * effectiveRadius ));
        speciesHandler.copyAndAddObject(species);
    }

    //sets up two particles with small overlap (such that forces are in equilibrium) at edge of domain, moving to the right
    void setupInitialConditions() override
    {
        auto s = speciesHandler.getLastObject();

        //add insertion boundaries
        SphericalParticle p;
        p.setSpecies(s);
        CubeInsertionBoundary b;
        b.setHandler(&boundaryHandler);
        b.set(p,1000,getMin(),{-1,1,1},{0,0,0},{0,0,0},.1,.1);
        b.setInitialVolume(1);
        boundaryHandler.copyAndAddObject(b);
    }
};

int main()
{
    InsertionBoundaryMPI2Test dpm;
    dpm.solve();
    return 0;
}
