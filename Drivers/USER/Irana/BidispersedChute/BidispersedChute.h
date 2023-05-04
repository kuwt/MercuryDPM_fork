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

#ifndef MERCURYDPM_BIDISPERSEDCHUTE_H
#define MERCURYDPM_BIDISPERSEDCHUTE_H

#include <iomanip>
#include <string.h>
#include <random>

#include "Chute.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "BidispersedChuteParameters.h"

class BidispersedChute : public Chute
{
public:
    
    BidispersedChute(const BidispersedChuteParameters& bidispersedChuteParameters = BidispersedChuteParameters());
    
    virtual ~BidispersedChute(){};
    
    //Since the time step is only written to the terminal when datafiles are written, we have to do it manually.
    void actionsBeforeTimeStep() override
    {
        //if (ceil(getTime()) != ceil(getTime() + getTimeStep()))
            //printTime();
    }
    
    //Set up periodic walls, rough bottom, add flow particles
    void setupInitialConditions() override;
    
    void setStandardDeviation(Mdouble sigma);
    
    //add flow particles
    void insertParticles(unsigned int numberOfLargeToGo, unsigned int numberOfSmallToGo);
    
    //defines type of flow particles
    void insertLargeParticle(unsigned int& numberOfLargeToGo);
    
    void insertSmallParticle(unsigned int& numberOfSmallToGo);
    
    void insertOneParticle(BaseParticle& p0);
    
    virtual void setBoundaries();
    
    void setChuteProperties();
    
    void setTimeAndSaveCount();
    
    virtual void setSpeciesProperties();
    
    void printTime() const override
    {
        logger(INFO, "t = %, tmax = %", getTime(), getTimeMax());
    }
    
    void makeLong()
    {
        isPeriodicInX = false;
    }
    
    void writeVTKFiles() const
    {
        DPMBase::writeVTKFiles();
    }
protected:
    BidispersedChuteParameters parameters;
    
    std::default_random_engine generator;
    std::normal_distribution<Mdouble> distribution;
    bool isPeriodicInX = true;
};

#endif //MERCURYDPM_BIDISPERSEDCHUTE_H
