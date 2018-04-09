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
#include "Species/SinterFrictionSpecies.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <sstream>
//#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
//#include "Walls/AxisymmetricIntersectionOfWalls.h"

/**
 * Test of the sintering properties of a single particle pair.
 * Two particles initially touching each other are sintered together.
 */
class Sinter : public Mercury3D {
public:


    void printTime() const override {
        static unsigned count = 0;
        if (++count==1) writeEneHeader(std::cout);
        writeEneTimestep(std::cout);
    }

    void writeEneHeader(std::ostream &os) const override {
        os << std::setw(13) << "time"
        << std::setw(13) << "overlap"
        << std::setw(13) << "plasticO"
        << std::setw(13) << "normalForce"
        << std::endl;
    }

    void writeEneTimestep(std::ostream &os) const override {
        auto i = dynamic_cast<const SinterInteraction *>(interactionHandler.getLastObject());
        if (i != nullptr) {
            os << std::setw(13) << getTime() << ' '
            << std::setw(13) << i->getOverlap() << ' '
            << std::setw(13) << i->getPlasticOverlap() << ' '
            << std::setw(13) << i->getForce().Z
            << std::endl;
        } else {
            logger(ERROR, "no SinterInteraction found");
        }
    }
};

int main(int argc UNUSED, char *argv[] UNUSED) {

    //create new simulation; as no specialised functions are needed,
    //Mercury3D is used directly instead of creating a derived class as usual.
    Sinter ps;
    Mdouble radius = 2e-6; //1 um
    Mdouble elasticModulus = 1e9;
    Mdouble rest = 0.5;

    //add Species
    auto s = ps.speciesHandler.copyAndAddObject(SinterFrictionSpecies());
    s->setDensity(1000);
    Mdouble mass = s->getMassFromRadius(radius);
    Mdouble meanOverlap = 0.1 * radius;
    Mdouble k = 4. / 3. * elasticModulus * std::sqrt(2.0 * radius * meanOverlap);
    logger(INFO, "k=%", k);
    s->setStiffnessAndRestitutionCoefficient(k, rest, mass);
    s->setPlasticParameters(k, 1.0001 * k, 0, 1);
    //s->setSinterAdhesion(2e-6/radius);//typical adhesion force 1e-6 N
    //s->setInverseSinterViscosity(1e-25);
    s->setSinterForceAndTime (2e-6/*N*/, 0.02/*s*/, radius);
    //s->setSinterRate(1e-9);
    Mdouble tc = s->getCollisionTime(mass);
    logger(INFO, "tc=%", tc);
    //3.88772*r/3.5e-1

    //add particle
    BaseParticle p;
    p.setSpecies(s);
    p.setRadius(radius);
    p.setPosition(p.getRadius() * Vec3D(0, 0, 1.0));
    ps.particleHandler.copyAndAddObject(p);

    //add second particle (with tiny overlap)
    p.setPosition(p.getRadius() * Vec3D(0, 0, -1.0 + 1e-13));
    ps.particleHandler.copyAndAddObject(p);

    //set time-stepping and output parameters
    ps.setFileType(FileType::ONE_FILE);
    ps.setXBallsAdditionalArguments(" -v0 -solidf ");
    ps.setGravity(Vec3D(0, 0, 0));
    ps.setTimeStep(0.02 * tc);
    ps.setSaveCount(100);
    ps.setTimeMax(1e4 * ps.getTimeStep());
    //ps.setSaveCount(1000);
    //ps.setTimeMax(100);
    ps.setMin(radius * Vec3D(-1, -1, -2));
    ps.setMax(radius * Vec3D(1, 1, 2));
    ps.setName("SinterParticlePairSelfTest");

    //run
    ps.solve();

    std::stringstream ss;
    ss << radius;
    helpers::writeToFile("SinterParticlePairSelfTest.gnu",
                         "set xlabel 't [s]'\n"
                          "r=" + ss.str() + "\n"
                          "p 'SinterParticlePairSelfTest.ene' u 1:($2/r) w lp t 'delta/r', '' u 1:($3/r) w lp t 'delta0/r'"
    );
    std::cout << "Execute 'gnuplot SinterParticlePairSelfTest.gnu' to view output" << std::endl;


}
