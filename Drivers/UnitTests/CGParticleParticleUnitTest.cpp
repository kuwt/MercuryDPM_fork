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

class TwoParticles : public Mercury3D
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
        p.setPosition({-0.05, 0, -0.95});
        particleHandler.copyAndAddObject(p);
        p.setPosition({0.05, 0, 0.95});
        particleHandler.copyAndAddObject(p);
    }
};

int main()
{
    logger(INFO," Simulates a particle-particle collision.\n"
            " Checks cg \n"
            "   - at one instant in time\n"
            "   - standard fields (density and stress)\n"
            "   - for O, Z, and XYZ coordinates\n"
            "   - for both Gauss and Lucy kernel functions\n");

    TwoParticles dpm;
    dpm.setName("CGParticleParticleUnitTest");
    dpm.setTimeStep(1e-4);
    dpm.setTimeMax(0);
    dpm.setFileType(FileType::NO_FILE);

    //Define coarse graining objects
    Mdouble width = 0.15;

    //Unresolved
    auto O = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());

    //Resolved in z at particle centre
    auto Z = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::Z>(width,1));
    Z->setMin({-1,-1,0.95});
    Z->setMax({1,1,0.95});

    //Resolved in z at particle centre with Gauss cg function
    auto G = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::Z,CGFunctions::Gauss>(width,1));
    G->setMin({-1,-1,0.95});
    G->setMax({1,1,0.95});

    //Resolved in z, at particle-particle contact
    auto C = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::Z>(width,1));

    //Resolved in xyz at particle centre
    auto XYZ = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::XYZ>(width,1));
    XYZ->setMin({-0.95,-1,0.95});
    XYZ->setMax({1.05,1,0.95});

    //Resolved in xz, at particle-particle contact
    auto cg2 = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::XZ,CGFunctions::Gauss>(3*width,30));
    auto cgG = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::XY,CGFunctions::Gauss>(width,50));
    auto cg = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::XY>(3*width,50));

    dpm.solve();

    //A few checks
    logger(INFO,"Checking a few cg parameters");

    //Checks density for unresolved cg, rho = M/V = 2/16 = 0.125
    helpers::check(O->getPoint(0).getDensity(),0.125,1e-15, "Average density");

    //Checks density at the center of one particle for z-resolved cg, rho = m_2*phi_2(0),
    //assuming the distance between the particles is bigger than the cutoff
    Mdouble phiL0 = 1.25 / width; //peak value of Lucy
    helpers::check(Z->getPoint(0).getDensity(),phiL0/4.0,1e-15, "Z-resolved density at particle centre");

    //Check the same for a Gauss cg function
    Mdouble phiG0 = 1.0 / (sqrt_2 * sqrt_pi * erf(3.0/sqrt_2))/width; //peak value of Gaussian cutoff at 3*width
    helpers::check(G->getPoint(0).getDensity(),phiG0/4.0,1e-15, "Z-resolved density for Gaussian CG");

    //Check the same for fully-resolved cg
    Mdouble phi30 = 105.0/16.0/pi / mathsFunc::cubic(width); //peak value of Lucy in 3D
    helpers::check(XYZ->getPoint(0).getDensity(),phi30,1e-11, "Fully-resolved density at particle centre");


    //Checks if the unresolved cg returns the right value for stress, sigma_zz = f_12*r_12/V
    Mdouble distance = Vec3D::getDistance(dpm.particleHandler.getObject(0)->getPosition(),dpm.particleHandler.getObject(1)->getPosition());
    Mdouble overlap = 2.0-distance;
    helpers::check(O->getPoint(0).getContactStress().trace(),2e5*overlap*distance/16,1e-11, "Average stress");

    //Checks if z-resolved cg returns the right value for stress at the center of the contact,
    //sigma_zz = f_12*r_12/r_12z/A, assuming the cutoff is smaller than the branch vector
    Mdouble distanceZ = dpm.particleHandler.getObject(1)->getPosition().Z-dpm.particleHandler.getObject(0)->getPosition().Z;
    helpers::check(C->getPoint(0).getContactStress().trace(),2e5*overlap*distance/distanceZ/4,1e-11, "Z-resolved stress at contact point");

    return 0;
}

