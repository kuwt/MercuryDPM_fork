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

#include<iostream>

#include "DPMBase.h"



class InsertionBoundaryTest : public DPMBase{
public:
	
	void setupInitialConditions()
	{
		wallHandler.clear();
		InfiniteWall w0;
		w0.set(Vec3D(-1, 0, 0), -getXMin());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 1, 0, 0),  getXMax());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0,-1, 0), -getYMin());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0, 1, 0),  getYMax());
		wallHandler.copyAndAddObject(w0);
        
        boundaryHandler.clear();
        InsertionBoundary b0;
        b0.set(Vec3D(getXMin(),getYMin(),getZMax()),Vec3D(getXMax(),getYMax(),getZMax()),Vec3D(-1.0,-1.0,0.0),Vec3D(1.0,1.0,0.0),0.05,0.15,100);
        boundaryHandler.copyAndAddObject(b0);
        
        particleHandler.clear();
        BaeParticle p0;
        p0.setPosition(Vec3D(0.5*(getXMin()+getXMax()),0.5*(getYMin()+getYMax()),0.5*(getZMin()+getZMax())));
        p0.setRadius(0.1);
        p0.set_Mass(1.0);                              
        p0.set_inertia(1.0);
        particleHandler.copyAndAddObject(p0);
	}
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
 	InsertionBoundaryTest problem;
 	problem.set_name("insertionBoundaryTest");
    problem.setGravity(Vec3D(0.0,0.0,0.0));
    problem.setXMax(1);
    problem.setYMax(1);
    problem.setTimeStep(1e-3);
    problem.speciesHandler.getObject(0)->setStiffness(1e5);
    problem.setDissipation(0);
    problem.setTimeMax(1);
    problem.setSaveCount(10);
    problem.setupInitialConditions();
    problem.boundaryHandler.get_Boundary(0)->checkBoundaryBeforeTimeStep(problem.particleHandler,problem.wallHandler,problem.random);
    problem.writeRestartFile();
    problem.write(std::cout,false);
    
    problem.solve();
}
