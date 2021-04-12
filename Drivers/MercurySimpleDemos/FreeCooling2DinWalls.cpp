
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

///todo{This code is not working as is wanted}
#include<iostream>
#include <Species/LinearViscoelasticSpecies.h>

#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include <ctime>
#include <Boundaries/PeriodicBoundary.h>
#include <Boundaries/CubeInsertionBoundary.h>
/// In this file 32^2 particles with the same velocity are placed in a bi-axial box. This makes them collide with the

/*#ifdef WITHGPERFTOOLS
    #include <gperftools/profiler.h>
#endif*/


class FreeCooling2DinWalls : public Mercury2D{
public:

    void setupInitialConditions() override {
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(10000);
        species->setDissipation(0.04);
        species->setStiffness(1e4);

        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(2e-4);
        p.setVelocity(Vec3D(1,1,1)*1e-3);

        CubeInsertionBoundary insertionBoundary;
        insertionBoundary.set(p,100,getMin(),getMax(),-p.getVelocity(),p.getVelocity(),p.getRadius(),p.getRadius());
        insertionBoundary.setInitialVolume(p.getVolume()*N);
        insertionBoundary.checkBoundaryBeforeTimeStep(this);
        setMeanVelocity({0,0,0});

        double mass  = species->getMassFromRadius(p.getRadius());
        double rest = species->getRestitutionCoefficient(mass);
        double tc = species->getCollisionTime(mass);
        logger(INFO,"Restitution %",rest);
        logger(INFO,"Collision time %",tc);
        logger(INFO,"N %",particleHandler.getNumberOfObjects());

        /*wallHandler.clear();
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D( 1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D( 0,-1, 0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D( 0, 1, 0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);*/

        PeriodicBoundary pb;
        pb.set(Vec3D(1,0,0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(pb);
        pb.set(Vec3D(0,1,0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(pb);
    }

    int N = 1;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{

	std::cout<<"In this file 32^2 particles with the same velocity are placed "
	"in a bi-axial box. This makes them collide with the walls and eachother. "
	"Afterwards the same run is performed with hgrid on. It tests the working "
	"(and speedup) of the hgrid."<<std::endl;

	FreeCooling2DinWalls problem;
	problem.setName("FreeCooling2DinWalls");

    problem.N=100;

    problem.setGravity(Vec3D(0.0,0.0,0.0));
    problem.setTimeStep(5e-5);
    problem.setSaveCount(2000);
    problem.setTimeMax(1000.0);
    problem.setMax({0.01,0.01,0.01});


    problem.setHGridMaxLevels(1);
	problem.setHGridCellOverSizeRatio(1.2);
	problem.setHGridUpdateEachTimeStep(false);
    problem.setFileType(FileType::ONE_FILE);

	problem.solve();
    return 0;
}
