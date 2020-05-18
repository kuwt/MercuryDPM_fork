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

#include <Walls/InfiniteWall.h>
#include "Mercury3D.h"
#include "CG/CG.h"
#include "Species/LinearViscoelasticSpecies.h"
using namespace constants;

class ParticleWall : public Mercury3D
{
public:
    void setupInitialConditions() override {
        //define a particle species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(0.75 / constants::pi); // such that mass = 1
        species->setStiffnessAndRestitutionCoefficient(2e5, 0.1, 1.0);

        //set domain size
        setDomain({-1,-1,-2},{1,1,2});

        //define two particles
        SphericalParticle p(species);
        p.setRadius(1.0);
        p.setPosition({0.05, 0, 0.95});
        particleHandler.copyAndAddObject(p);

        InfiniteWall w({-0.05, 0, -0.95},{0,0,0},species);
        wallHandler.copyAndAddObject(w);
    }
};

int main()
{
    logger(INFO," Simulates a particle-wall collision.\n"
            " Checks cg \n"
            "   - at one instant in time\n"
            "   - standard fields (density and stress)\n"
            "   - for O, Z, and XYZ coordinates\n"
            "   - for both Gauss and Lucy kernel functions\n");

    ParticleWall dpm;
    dpm.setName("CGParticleWallUnitTest");
    dpm.setTimeStep(1e-4);
    dpm.setTimeMax(0);
    dpm.setFileType(FileType::NO_FILE);

    //Define coarse graining objects
    Mdouble width = 0.075;

    //Unresolved
    auto O = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());

    //Resolved in z, at particle center
    auto Z = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::Z>(width,1));
    Z->setMin({-1,-1,0.95});
    Z->setMax({1,1,0.95});

    //Resolved in z with Gauss cg function
    auto G = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::Z,CGFunctions::Gauss>(width,1));
    G->setMin({-1,-1,0.95});
    G->setMax({1,1,0.95});

    //Resolved in z, on contact line
    auto C = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::Z>(width,1));
    C->setMin({-1,-1,0.25});
    C->setMax({1,1,0.25});

    //Resolved in xyz, at particle center
    auto XYZ = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::XYZ>(width,1));
    XYZ->setMin({-0.95,-1,0.95});
    XYZ->setMax({1.05,1,0.95});

    dpm.solve();

    //A few checks
    logger(INFO,"Checking a few cg parameters");

    //Checks if the unresolved cg returns the right value for stress, sigma_zz = f_12*r_12/V
    Vec3D contactPoint = dpm.interactionHandler.getObject(0)->getContactPoint();
    Mdouble overlap = dpm.interactionHandler.getObject(0)->getOverlap();
    Mdouble distance = dpm.interactionHandler.getObject(0)->getDistance();
    Mdouble distance2 = Vec3D::getDistance(contactPoint,dpm.particleHandler.getObject(0)->getPosition());
    Mdouble distance3 = distance+0.5*overlap;
    logger(INFO,"Distance % % % Overlap % Contact % Distance %",distance,distance2,distance3,overlap,contactPoint);
    helpers::check(O->getPoint(0).getContactStress().trace(),2e5*overlap*distance2/16,1e-11, "Average stress");
    ///\todo TW should there be a discrepancy between distance and distance2?

    //Checks if z-resolved cg returns the right value for stress at the center of the contact,
    //sigma_zz = f_12*r_12/r_12z/A, assuming the cutoff is smaller than the branch vector
    Mdouble distanceZ = dpm.particleHandler.getObject(0)->getPosition().Z;
    helpers::check(C->getPoint(0).getContactStress().trace(),2e5*overlap*distance/distanceZ/4,1e-11, "Z-resolved stress at contact point");


    return 0;
}

