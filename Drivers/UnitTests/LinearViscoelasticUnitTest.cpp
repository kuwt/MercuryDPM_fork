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

#include <Mercury3D.h>

//! [T11:contactModel]
#include <Species/LinearViscoelasticSpecies.h>
//! [T11:contactModel]

///This code tests the linear viscoelastic  behavior of two particles.
class viscoElasticUnitTest : public DPMBase {
public:

    explicit viscoElasticUnitTest(Mdouble radius) {

        //------------------
        //Global parameters
        setName("LinearViscoElasticUnitTest");
        setFileType(FileType::ONE_FILE);
        setGravity(Vec3D(0.0, 0.0,0.0));

        //-------------------
        //Boundary parameters
        setXMax(2.0 * radius);
        setYMax(radius);
        setZMax(radius);
        setXMin(-getXMax());
        setYMin(-getYMax());
        setZMin(-getZMax());

        //-------------------
        //Setup parameters:
        const Mdouble density = 1000;
        const Mdouble k1 = 1.0 * radius; //[Stiffness depends on particle radius]
        const Mdouble restitutionCoefficient = 1.0;
        const Mdouble pDepth = 0.1;
        //-------------------

        //------------------
        //Species:
        LinearViscoelasticSpecies sf;
        sf.setHandler(&speciesHandler);
        sf.setDensity(density);

        const Mdouble mass = sf.getMassFromRadius(radius);
        sf.setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, mass);

        const Mdouble collisionTime = sf.getCollisionTime(mass);
        logger(INFO, "Collision time %", collisionTime);

        //-------------------
        auto species = speciesHandler.copyAndAddObject(sf);

        //-------------------
        //Particle properties:
        SphericalParticle P0, P1;
        //SphericalParticle P0, P1;
        P0.setSpecies(species);
        P1.setSpecies(species);

        P0.setRadius(radius);
        P1.setRadius(radius);

        P0.setPosition(Vec3D(-1.3 * radius, 0, 0));
        P0.setVelocity(Vec3D(0.001,0,0));
        P1.setPosition(-P0.getPosition());
        P1.setVelocity(-P0.getVelocity());

        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);

        //------------------
        // Time
        setSaveCount(10000);
        setTimeMax(10);
        setTimeStep(0.0000005 * 50);
        //------------------
    }

    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    viscoElasticUnitTest sp0(0.01);
    sp0.solve();

    //This helper is to see the Fn vs Overlap, rate of overlap.
    std::cout << "Execute 'load LinearViscoElasticUnitTest.gnu' to view output" << std::endl;
    helpers::writeToFile("LinearViscoElasticUnitTest.gnu",
                         "set xlabel 'displacement [{/Symbol d}]'\n"
                         "set ylabel 'force [f^n]'\n"
                         "set grid\n"
                         "plot 'LinearViscoElasticUnitTest.fstat' u 7:9 w lp"
    );

    return 0;

}
