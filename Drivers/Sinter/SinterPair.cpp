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

//! [St:headers]
#include "Mercury3D.h"
#include "Species/SinterSpecies.h"
//! [St:headers]

using std::cout;
using std::endl;

/// Single particle pair, sintered slowly.

//! [St:class]
class SinterPair : public Mercury3D
{
public:
    explicit SinterPair (Mdouble radius)
    {
        std::string r = helpers::to_string(radius);
        setName("SinterPair"+r);
        helpers::writeToFile("SinterPair"+r+".gnu",
                             "set xlabel 'time [s]'\n"
                              "set ylabel 'x/a [nm]'\n"
                              "r="+r+"\n"
                              "plot 'SinterPair"+r+".fstat' u ($1):(sqrt($7/r)) w lp\n"
        );

        setGravity({0,0,0});
        setSaveCount(20000);
        setTimeMax(20);
        setMin({-2*radius,-radius,-radius});
        setMax({2*radius,radius,radius});
        setParticlesWriteVTK(true);

        const Mdouble stiffness = 1e-2 * radius;
        const Mdouble restitutionCoefficient = 0.1;
        const Mdouble density = 1005;

        //! [St:speciesProp]
        SinterSpecies s;
        s.setHandler(&speciesHandler);
        s.setDensity(density);
        const Mdouble mass = s.getMassFromRadius(radius);
        s.setStiffnessAndRestitutionCoefficient(stiffness,restitutionCoefficient,mass);
        s.setPlasticParameters(stiffness, 10*stiffness, stiffness, 0.16);
        const Mdouble collisionTime = s.getCollisionTime(mass);
        logger(INFO,"Collision time %",collisionTime);
        s.setSinterType(SINTERTYPE::CONSTANT_RATE);
        s.setSinterRate(4e-10/radius);
        s.setSinterAdhesion(0.013*stiffness);
        //adhesiveForce = sinterAdhesion*radius;
        auto species = speciesHandler.copyAndAddObject(s);
        //! [St:speciesProp]

        setTimeStep(0.02*collisionTime);

        //! [St:createParticle]
        SphericalParticle p;
        p.setSpecies(species);
        p.setRadius(radius);
        p.setPosition({-(1-1e-15)*radius,0,0});
        particleHandler.copyAndAddObject(p);
        p.setPosition(-p.getPosition());
        particleHandler.copyAndAddObject(p);
        //! [St:createParticle]
    }

    void printTime() const override
    {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
         << ", r=" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getLastObject()->getRadius()
         //<< ", Ekin=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
                  << std::endl;
        std::cout.flush();
    }

};
//! [St:class]

//! [St:main]
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //! [St:problemSetup]
    SinterPair sp0(2e-6);
    sp0.solve();
    SinterPair sp1(1.5e-6);
    sp1.solve();
    SinterPair sp2(5e-7);
    sp2.solve();
    //! [St:problemSetup]

    //! [St:output]
    helpers::writeToFile("SinterPair.gnu",
                         "set xlabel 'time [s]'\n"
                          "set ylabel 'x/a'\n"
                          "plot [0:200] 'SinterPair5e-07.fstat' u ($1):(sqrt($7/5e-07)) w lp lt rgb 'royalblue'\n"
                          "replot 'SinterPair1.5e-06.fstat' u ($1):(sqrt($7/1.5e-06)) w lp lt rgb 'light-red'\n"
                          "replot 'SinterPair2e-06.fstat' u ($1):(sqrt($7/2e-06)) w lp lt rgb 'sea-green'"
    );
    std::cout << "Execute 'gnuplot SinterPair.gnu' to view output" << std::endl;
    //! [St:output]
}
//! [St:main]
