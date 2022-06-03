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

//! [Lees:headers]
#include "Mercury2D.h"
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
//! [Lees:headers]

//! [Lees:class]
class LeesEdwardsSelfTest : public Mercury2D
{
    public:

    void setupInitialConditions() override {
        //set parameters to define the species properties

        setName("LeesEdwardsSelfTest");
        dataFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::ONE_FILE);

        const Mdouble particleRadius = 0.5; // such that diameter is one
        const Mdouble collisionTime = 0.05; //relatively stiff particles
        const Mdouble restitution = 0.95; //restitution close to glass particles
        const Mdouble sizeDistribution = 1.5;
        const Mdouble volumeFraction = 0.82;
        const Mdouble velocity = 15.0;

        //define species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(4.0 / constants::pi); //such that particle mass is one
        Mdouble mass = species->getMassFromRadius(particleRadius);
        species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(
                collisionTime, restitution, restitution, mass); //set stiffness and dissipation
        species->setSlidingFrictionCoefficient(0.5);

        setTimeStep(0.02 * species->getCollisionTime(mass));
        setTimeMax(4.0); //run until the situation is static
        setSaveCount(100);

        //set domain size
        setMin({0, 0, 0});
        setMax({15, 15, 1});
        //set gravity
        setGravity({0, 0, 0});


        //define leesEdwardsBoundary
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.set(
                [velocity](double time) { return time * velocity; },
                [velocity](double time UNUSED) { return velocity; },
                getXMin(), getXMax(), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(leesEdwardsBoundary);

        //define common particle properties
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(particleRadius);
        Mdouble rMin = 2.0 * particleRadius / (sizeDistribution + 1);
        Mdouble rMax = sizeDistribution * rMin;
        p.setRadius(rMax);
        logger(INFO, "Inserting particles of diameter % to %, volumeFraction %", 2.0 * rMin, 2.0 * rMax,
               volumeFraction);

        Mdouble particleVolumeMax = volumeFraction * (getXMax() - getXMin()) * (getYMax() - getYMin());
        Mdouble particleVolume = 0;
        Vec3D position = {0, 0, 0};
        while (particleVolume + 0.5 * p.getVolume() < particleVolumeMax) {
            position.X = random.getRandomNumber(getXMin(), getXMax());
            position.Y = random.getRandomNumber(getYMin(), getYMax());
            p.setPosition(position);
            particleHandler.copyAndAddObject(p);
            particleVolume += p.getVolume();
            p.setRadius(random.getRandomNumber(rMin, rMax));
        }

    }


};
//! [Lees:class]

int main()
{
    //instantiate the class
    LeesEdwardsSelfTest problem;
    //set output and time stepping properties
    problem.setXBallsAdditionalArguments("-w0 -v0 -solidf -cmode 5");
    //solve
    problem.solve();

    logger(INFO,"Execute 'source LeesEdwardsSelfTest.sh' to get coarse-grained statistics");
    helpers::writeToFile("LeesEdwardsSelfTest.sh","../MercuryCG/fstatistics LeesEdwardsSelfTest -stattype XY -w 0.25 -h 0.25 -tmin 2\n"
     "../MercuryCG/fstatistics LeesEdwardsSelfTest -stattype Y -w 0.25 -h 0.1 -tmin 2 -o LeesEdwardsSelfTest.Y.stat");

    logger(INFO,"Run 'LeesEdwardsSelfTest.m' in MATLAB/octave to visualise the statistical output");
    helpers::writeToFile("LeesEdwardsSelfTest.m","%% 2D velocity field v_x(x,y)\n"
     "addpath('../MercuryCG/')\n"
     "data = loadstatistics('LeesEdwardsSelfTest.stat');\n"
     "colormap(jet)\n"
     "contourf(data.x,data.y,data.VelocityX,20,'EdgeColor','none')\n"
     "c = colorbar\n"
     "c.Label.String = '\\rho';\n"
     "title('Velocity')\n"
     "xlabel('x')\n"
     "ylabel('z');\n"
     "axis equal\n"
     "%% 1D velocity field v_x(y)\n"
     "dataY = loadstatistics('LeesEdwardsSelfTest.Y.stat');\n"
     "plot(dataY.y,dataY.VelocityX)\n"
     "xlabel('y')\n"
     "ylabel('v_x');\n"
     "axis equal");

    return 0;
}
//! [Lees:main]
