//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
#include "DPMBase.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include <Logger.h>

/*!
 * \brief Compresses 2 particles (vertically) until they have an extreme overlap 
 * \details
This code tries to emulate the following inputfiles from Stefan
c3d.ini:
 2 0 0 0 0 .01 .01 .01
.005 .005 .009 0. 0 0 .0010 0
.005 .005 .004 0. 0 0 .0040 0

and par.ini:
0 1 0
1e-2 1e-7
.0001 .0002
2000 1e5 1e-2 0
0 0

and walls.ini
 1 10 0 0.01
 2 10 0 0.01
 3 10 0 0.01
 4 10 0 0.01
 5 10 0 0.01
 6 10 .005 0.02
 */
class ExtremeOverlapWithWallsUnitTest : public DPMBase{
public:

    void setupInitialConditions() override {
        Mdouble r0=0.004, r1=0.001;

        setXMax(2.0*(r0+r1));
        setXMin(0.0);
        setYMax(r0);
        setYMin(-r0);
        setZMax(r0);
        setZMin(-r0);

        ///\bug setting gravity here was a quick fix when I changed the default gravity to zero. Should gravity be 0 here?
        setGravity(Vec3D(0,0,0));

        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1,0,0), Vec3D(0,0,0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1,0,0),  Vec3D(getXMax(),0,0));
        //should be set to 0.005*sin(2*pi*t/0.02))
        //w0.setVelocity(Vec3D(0.0,0.0,-constants::pi/2.0));
		w0.setPrescribedPosition([this,r1] (double time)
        {
            return Vec3D(getXMax() - 4.0*r1 * std::sin(time * 2.0 * constants::pi/0.02),0,0);
        });
        wallHandler.copyAndAddObject(w0);

        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setPosition(Vec3D(r0,0,0));
        p0.setVelocity(Vec3D(0.0,0.0,0.0));
        p0.setRadius(r0);
        particleHandler.copyAndAddObject(p0);
        p0.setPosition(Vec3D(r0+r0+r1,0,0));
        p0.setVelocity(Vec3D(0.0,0.0,0.0));
        p0.setRadius(r1);
        particleHandler.copyAndAddObject(p0);
    }

//    void printTime() const override
//	{
//        if (interactionHandler.getNumberOfObjects()>1)
//        std::cout
//        << getTime() << " "
//        << particleHandler.getObject(0)->getPosition().Z << " "
//        << particleHandler.getObject(1)->getPosition().Z << " "
//        << interactionHandler.getObject(1)->getContactPoint().Z <<std::endl;
//    }

};

int main(int argc, char *argv[])
{
	///Start off my solving the default problem
	ExtremeOverlapWithWallsUnitTest OverlapProblem;
    auto species = OverlapProblem.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    species->setDensity(2000);
    species->setStiffness(1e5);
    species->setDissipation(1e-2);
    OverlapProblem.setName("ExtremeOverlapWithWallsUnitTest");
    OverlapProblem.setFileType(FileType::NO_FILE);
	OverlapProblem.dataFile.setFileType(FileType::ONE_FILE);
    //OverlapProblem.setWallsWriteVTK(true);
    //OverlapProblem.setParticlesWriteVTK(true);
	OverlapProblem.setSaveCount(1000);
    OverlapProblem.setDimension(3);
    OverlapProblem.setTimeStep(1e-7);
	OverlapProblem.setTimeMax(0.01);
    // OverlapProblem.dataFile.getFstream().precision(10);
	OverlapProblem.solve(argc,argv);
    
    Vec3D position = OverlapProblem.particleHandler.getObject(0)->getPosition();
    helpers::check(position,Vec3D(0.00402345, 0, 0), 1e-8,"First particle position");

    position = OverlapProblem.particleHandler.getObject(1)->getPosition();
    helpers::check(position,Vec3D(0.00901182, 0, 0), 1e-8,"Second particle position");

    return 0;
}
