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


#include "BidispersedBelt.h"

void BidispersedBelt::setupInitialConditions()
{
    BidispersedChute::setupInitialConditions();
    prescribeBeltParticleVelocity();
}

void BidispersedBelt::prescribeBeltParticleVelocity()
{
    for (BaseParticle* p : particleHandler)
    {
        if (p->isFixed())
        {
            p->setPrescribedVelocity([this](Mdouble time){return Vec3D(this->getBeltSpeed(), 0, 0);});
        }
    }
}

Mdouble BidispersedBelt::getBeltSpeed()
{
    return beltSpeed_;
}

void BidispersedBelt::setBeltSpeed(Mdouble beltSpeed)
{
    beltSpeed_ = beltSpeed;
}

void BidispersedBelt::setBoundaries()
{
    BidispersedChute::setBoundaries();
    if (!isPeriodicInX_) //add walls on both sides of the domain
    {
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        //wall on left side of the domain
        w0.set(Vec3D(-1.0, 0.0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        //wall on the right side of the domain
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
    }
}

void BidispersedBelt::setPeriodicInX(bool isPeriodicInX)
{
    isPeriodicInX_ = isPeriodicInX;
}

void BidispersedBelt::setAutomaticBeltSpeed(bool isAutomaticSpeed)
{
    isAutomaticSpeed_ = isAutomaticSpeed;
    if (isPeriodicInX_)
    {
        logger(WARN, "cannot regulate belt-speed when using periodic boundary conditions in x-direction");
        isAutomaticSpeed_= false;
    }
}

///Adapt the speed of the belt, copied from Kasper's Belt
///Probably: keep COM in x-direction in the middle of the domain.
void BidispersedBelt::actionsBeforeTimeStep()
{
    BidispersedChute::actionsBeforeTimeStep();
    if (isAutomaticSpeed_)
    {
        static unsigned int beltSpeedUpdateTimer = 0;
        beltSpeedUpdateTimer++;
    
        if (beltSpeedUpdateTimer > 50000)
        {
        
            beltSpeedUpdateTimer = 0;
        
            double S = 0;
        
            for (BaseParticle* p : particleHandler)
            {
                if (p->isFixed())
                {
                    continue;
                }
            
                if (p->getIndSpecies() == 2)
                {
                    S = S + p->getPosition().X * pow(parameters.getSmallParticleRadius(), 3);
                }
                else if (p->getIndSpecies() == 1)
                {
                    S = S + p->getPosition().X * pow(parameters.getLargeParticleRadius(), 3);
                }
            }
        
        
            S = S / ((parameters.getNumberOfSmallParticles(getChuteLength() * getChuteWidth()) * pow(parameters.getSmallParticleRadius(), 3))
                     + (parameters.getNumberOfLargeParticles(getChuteWidth() * getChuteLength()) * pow(parameters.getLargeParticleRadius(), 3)));
        
            double dBelt = -(S - getXMax() / 2) / 2000;
        
            setBeltSpeed(getBeltSpeed() + dBelt);
            std::cout << "BS = " << getBeltSpeed() << std::endl;
            std::cout << "dBelt = " << dBelt << std::endl;
            std::cout << "COM = " << S << std::endl;
            prescribeBeltParticleVelocity();
        }
    }
    
}
