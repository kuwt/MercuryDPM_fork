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

// This code generates the plot in boundary statistics paper

//! [FP:headers]
#include <Mercury3D.h>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
//! [FP:headers]

//! [FP:class]
class FiveParticles : public Mercury3D
{
public:
	void setupInitialConditions() override {
		//set parameters to define the species properties
		const Mdouble particleRadius=0.5; // such that diameter is one
		const Mdouble collisionTime = 0.005; //relatively stiff particles
		const Mdouble restitution = 0.88; //restitution close to glass particles

		//define species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
		species->setDensity(6.0/constants::pi); //such that particle mass is one
		species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(collisionTime, restitution, restitution, species->getMassFromRadius(particleRadius)); //set stiffness and dissipation
		species->setSlidingFrictionCoefficient(0.5);

		//set gravity, time step
		setGravity(Vec3D(0,0,-1)); // such that gravity is one
		setTimeStep(0.02 * collisionTime);

		//set domain size
		setXMin(0); setYMin(0); setZMin(0);
		setXMax(5); setYMax(1); setZMax(2.5);

		//define common particle properties
		SphericalParticle p0;
		p0.setSpecies(speciesHandler.getObject(0));
		p0.setRadius(particleRadius);

        //define five fixed particles
		p0.fixParticle();
		p0.setPosition(Vec3D(0.5,0.5,0.1));
		particleHandler.copyAndAddObject(p0);
		p0.setPosition(Vec3D(1.5,0.5,-.2));
		particleHandler.copyAndAddObject(p0);
		p0.setPosition(Vec3D(2.5,0.5,0.0));
		particleHandler.copyAndAddObject(p0);
		p0.setPosition(Vec3D(3.5,0.5,0.1));
		particleHandler.copyAndAddObject(p0);
		p0.setPosition(Vec3D(4.5,0.5,0.05));
		particleHandler.copyAndAddObject(p0);

        //define five free particles
		p0.unfix();
		p0.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0.0,random.getRandomNumber(-.1,.1)));
		p0.setPosition(Vec3D(1.,0.5,1.0));
		particleHandler.copyAndAddObject(p0);
		p0.setPosition(Vec3D(2.0,0.5,1.0));
		particleHandler.copyAndAddObject(p0);
		p0.setPosition(Vec3D(3.0,0.5,1.0));
		particleHandler.copyAndAddObject(p0);
		p0.setPosition(Vec3D(4.0,0.5,1.0));
		particleHandler.copyAndAddObject(p0);
		p0.setPosition(Vec3D(1.5,0.5,2.0));
		particleHandler.copyAndAddObject(p0);
	}
};
//! [FP:class]

//! [FP:main]
int main(int argc, char *argv[])
{
	//instantiate the class
	FiveParticles fiveParticles;
	//set name
	fiveParticles.setName("FiveParticles");

	//set output and time stepping properties
	fiveParticles.setTimeMax(20); //run until the situation is static
	fiveParticles.setSaveCount(std::numeric_limits<unsigned >::max()); //save only the first and last time step
	//solve
	fiveParticles.solve(argc,argv);

    logger(INFO,"Execute 'source FiveParticles.sh' to get coarse-grained statistics of the last time step");
    helpers::writeToFile("FiveParticles.sh","../MercuryCG/fstatistics FiveParticles -stattype XZ -w 0.1 -h 0.05 -tmin 1 -tmax 30");

    logger(INFO,"Run 'FiveParticles.m' in MATLAB/octave to visualise the statistical output");
    helpers::writeToFile("FiveParticles.m","addpath('../../MercuryCG/')\n"
     "data = loadstatistics('FiveParticles.stat');\n"
     "colormap(1-gray)\n"
     "contourf(data.x,data.z,data.Density,20,'EdgeColor','none')\n"
     "c = colorbar\n"
     "c.Label.String = '\\rho';\n"
     "title('Density')\n"
     "xlabel('x')\n"
     "ylabel('z');\n"
     "axis equal\n"
     "%%\n"
     "particles=importdata('FiveParticles.data',' ',12);\n"
     "x=particles.data(:,1);\n"
     "z=particles.data(:,3);\n"
     "r=particles.data(:,7);\n"
     "a=linspace(0,2*pi,40);\n"
     "xCircle = sin(a);\n"
     "zCircle = cos(a);\n"
     "hold on;\n"
     "for i=1:length(x)\n"
     "  plot(x(i)+r(i)*xCircle,z(i)+r(i)*zCircle,'Color',.8*[1 1 1])\n"
     "end\n"
     "hold off");
}
//! [FP:main]
