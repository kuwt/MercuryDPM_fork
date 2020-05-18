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

#include<iostream>
#include <Species/LinearViscoelasticSpecies.h>

#include "Mercury2D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"

class my_problem : public Mercury2D
{
public:

    void setupInitialConditions()
    {
        //Put the particles on a grid with small random velocities
        double x = 0.6;
        double y = 0.6;
        double z = getZMin() + 2.;
        for (int i = 0; i < N; i++)
        {
            x += 1.2;
            if (x > getXMax() - 0.5)
            {
                x = 0.6;
                y += 1.2;
            }
            if (y > getYMax() - 0.5)
            {
                x = 0.6;
                y = 0.6;
                z += 1.2;
            }
            particleHandler.getObject(i)->setPosition(Vec3D(x, y, z));
            particleHandler.getObject(i)->setVelocity(Vec3D(drand48() * 0.01, drand48() * 0.01, drand48() * 0.01));
            particleHandler.getObject(i)->setRadius(0.5);

        }

        //One wall at the bottom and periodic in other directions
        InfiniteWall w0;
        w0.set(Vec3D(0.0, 0.0, -1.0), getMin());
        w0.setPrescribedPosition([this] (double time)
        {
            double t = getTime() - 1.0;
            if (t > 0.0)
            {
                return Vec3D(0.0,0.0,getZMin() + 2.0 * std::sin(t * 0.5));
            }
            else
            {
                return Vec3D(0.0,0.0,getZMin());
            }
        });
        wallHandler.copyAndAddObject(w0);
        PeriodicBoundary b0;
        b0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0, 1, 0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);

        //Set gravity to be down magitude 1
        setGravity(Vec3D(0.0, 0.0, -1.0));
    }

    unsigned int N;

///todo{DK: This is the old moving wall functionallity, new one has to be tested}
/*
protected:
void actionsBeforeTimeStep()
{
//After t=1.0 start to move the bottom wall
double t=getTime()-1.0;
if (t>0.0)
	{
	  wallHandler.getObject(0)->move(getZMin()-2.0*sin(t*0.5));
	}

}
*/
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc UNUSED, char *argv[] UNUSED)
{

    ///Start off my solving the default problem
    my_problem problem;
    problem.setSystemDimensions(3);
    problem.N = 2500;
    problem.setName("leidenfrost");
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    species->setDissipation(0.01);
    species->setStiffness(1e3);
    problem.setTimeStep(1e-3);
    problem.setSaveCount(50);
    species->setDensity(1);
    problem.setTimeMax(200);
    problem.setXMax(50);
    problem.setYMax(5);
    problem.setZMax(100);
    problem.solve();

}
