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
#include <vector>
#include <Walls/InfiniteWall.h>
#include "DPMBase.h"
#include "Mercury3D.h"
#include "Chute.h"

///This code examines the flow front of rough-bottom chute flow. The 
///flow is initialised on \f$x\in[0,FlowLength]\f$. Then the chute gradually 
///expands to the right as the flow develops, and is cut on the left 
///to minimize computation.
class FlowFrontChute : public Chute
{
public:
    void stretch()
    {
        ///prolong the chute 10-fold
        int stretchFactorBottom = 4;
        int stretchFactorFlow = 2;
        
        unsigned numberOfParticles = particleHandler.getSize();
        //add new particles
        for (unsigned i = 0; i < numberOfParticles; ++i)
        {
            BaseParticle* p  = particleHandler.getObject(i);
            if (!p->isFixed())
            {
                for (int j = 0; j < stretchFactorFlow; j++) //add flow particles
                {
                    Vec3D newPosition = p->getPosition() + Vec3D((getXMax() - getXMin()), 0, 0);
                    p->setPosition(newPosition);
                    particleHandler.copyAndAddObject(p);
                }
            }
            else
            {
                for (int j = 0; j < stretchFactorBottom; j++) //add bottom particles
                {
                    Vec3D newPosition = p->getPosition() + Vec3D((getXMax() - getXMin()), 0, 0);
                    p->setPosition(newPosition);
                    particleHandler.copyAndAddObject(p);
                }
            }
        }
        
        //stretch domain
        setXMax(getXMin() + stretchFactorBottom * (getXMax() - getXMin()));
        for (BaseWall* w : wallHandler)
        {
            InfiniteWall* w0 = dynamic_cast<InfiniteWall*> (w);
            if (w0 != nullptr && w0->getNormal().X == 1.0)
                w0->setPosition(Vec3D(getXMax(), 0.0, 0.0));
        }
        
    }
    
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    
    FlowFrontChute problem;
    problem.setName("FlowFrontChute");
    problem.stretch();
    
    problem.write(std::cout, false);
    problem.solve();
    
    problem.write(std::cout, false);
    problem.writeRestartFile();
}
