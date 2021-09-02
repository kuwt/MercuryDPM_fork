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

#include "PolydisperseInsertionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Math/RNG.h"
//#include <cassert>

/*!
 * \details Deprecated boundary which was used to insert PSDs into Mercury.
 * \deprecated Should be gone by Mercury 2.0. Instead, use the PSD class.
 */
/* Constructor */
PolydisperseInsertionBoundary::PolydisperseInsertionBoundary()
{
    logger(INFO, "In PolydisperseInsertionBoundary constructor");
    /* std::vector does all the memory allocation for you. 
     * By default, these are vectors of zero length. */
    logger(INFO, "About to leave PolydisperseInsertionBoundary copy constructor");
}

/* Copy constructor */
PolydisperseInsertionBoundary::PolydisperseInsertionBoundary(const PolydisperseInsertionBoundary& other)
        : InsertionBoundary(other)
{
    /* The new PolydisperseInsertionBoundary's generanda_ vector should point to
     * its own copies of each of the generanda_, so we need to copy those
     * across. */
    for (int i = 0; i < generanda_.size(); i++)
        generanda_[i] = other.generanda_[i]->copy();
    
    posMin_ = other.posMin_;
    posMax_ = other.posMax_;
    velMin_ = other.velMin_;
    velMax_ = other.velMax_;
    probabilitates_ = other.probabilitates_;
    sizeDispersities_ = other.sizeDispersities_;
    numbersInserted_ = other.numbersInserted_;
    massesInserted_ = other.massesInserted_;
    volumesInserted_ = other.volumesInserted_;
}

/* Destructor */
PolydisperseInsertionBoundary::~PolydisperseInsertionBoundary()
{
    // JMFT: Do we need to delete the elements of generanda_?
    for (auto p : generanda_)
        
        delete p;
}

PolydisperseInsertionBoundary* PolydisperseInsertionBoundary::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    logger(INFO, "PolydisperseInsertionBoundary::copy() const finished");
#endif
    return new PolydisperseInsertionBoundary(*this);
}

void PolydisperseInsertionBoundary::setGeometry(int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax)
{
    setMaxFailed(maxFailed);
    posMin_ = posMin;
    posMax_ = posMax;
    velMin_ = velMin;
    velMax_ = velMax;
}

void PolydisperseInsertionBoundary::addGenerandum(BaseParticle* generandum, double probability, double sizeDispersity)
{
    // Give the PolydisperseInsertionBoundary its own copy of the generandum.
    generanda_.push_back(generandum->copy());
    probabilitates_.push_back(probability);
    sizeDispersities_.push_back(sizeDispersity);
    numbersInserted_.push_back(0);
    massesInserted_.push_back(0);
    volumesInserted_.push_back(0);
    logger(INFO, "PolydisperseInsertionBoundary: added a new generandum, now have %. New generandum has weighting %",
           generanda_.size(), probability);
}


void PolydisperseInsertionBoundary::setGenerandum(unsigned int spec, BaseParticle* generandum, double probability,
                                                  double sizeDispersity)
{
    if (spec < generanda_.size())
        logger(INFO, "Setting the %-th species of a PolydisperseInsertionBoundary that so far has % species",
               spec, generanda_.size());
    else
        logger(ERROR,
               "Setting the %-th species of a PolydiserseInsertionBoundary with only % species is illegal. Use addGenerandum instead.",
               spec, generanda_.size());
    
    if (probability == 0)
        logger(WARN, "PolydisperseInsertionBoundary: Are you sure you want to set the probability to be 0?");
    
    generanda_[spec] = generandum->copy();
    probabilitates_[spec] = probability;
    sizeDispersities_[spec] = sizeDispersity;
    
    // Reset the counters for number, mass and volume
    numbersInserted_[spec] = 0;
    massesInserted_[spec] = 0;
    volumesInserted_[spec] = 0;
}

BaseParticle* PolydisperseInsertionBoundary::generateParticle(RNG& random)
{
    /* First choose what particle species to generate. */
    Mdouble totalprob = 0;
    for (auto p : probabilitates_)
        totalprob += p;
    
    Mdouble check = random.getRandomNumber(0, totalprob);
    unsigned int spec;
    logger(VERBOSE, "PolydisperseInsertionBoundary: check = % out of %",
           check, totalprob);
    for (int i = 0; i < generanda_.size(); i++)
    {
        if (check >= probabilitates_[i])
            check -= probabilitates_[i];
        else
        {
            spec = i;
            break;
        }
    }
    
    auto P = generanda_[spec]->copy();
    
    /* The 'reference' particle for this species has a radius, but we allow some
     * sizeDispersity. */
    double radius = P->getRadius()
                    * random.getRandomNumber(1 - sizeDispersities_[spec], 1 + sizeDispersities_[spec]);
    P->setRadius(radius);
    
    /* JMFT: TODO: These do *not* give the correct values! 
     * They give the number &c. of each species that the
     * PolydisperseInsertionBoundary has _attempted_ to place, not the number
     * that have actually been placed successfully. */
    numbersInserted_[spec]++;
    massesInserted_[spec] += P->getMass();
    // volumesInserted_[spec] += P->getVolume();
    
    return P;
}

/* JMFT: TODO: We should think how to recycle this code from
 * CubeInsertionBoundary more efficiently */
void PolydisperseInsertionBoundary::placeParticle(BaseParticle* P, RNG& random)
{
    Vec3D pos, vel;
    pos.X = random.getRandomNumber(posMin_.X, posMax_.X);
    pos.Y = random.getRandomNumber(posMin_.Y, posMax_.Y);
    pos.Z = random.getRandomNumber(posMin_.Z, posMax_.Z);
    vel.X = random.getRandomNumber(velMin_.X, velMax_.X);
    vel.Y = random.getRandomNumber(velMin_.Y, velMax_.Y);
    vel.Z = random.getRandomNumber(velMin_.Z, velMax_.Z);
    P->setPosition(pos);
    P->setVelocity(vel);
}

// JMFT: TODO
void PolydisperseInsertionBoundary::read(std::istream& is)
{
    InsertionBoundary::read(is);
    /*
    for (int i = 0; i < generanda_.size(); i++)
    {
        BaseParticle* particleToCopy = new SphericalParticle;
        // BaseParticle::write writes the extra word 'BaseParticle', which will be
        // ignored by BaseParticle::read. To avoid an off-by-one error, we need to
        // get rid of this extra word first...
        std::string dummy;
        is >> dummy;

        // Now read the particleToCopy.
        particleToCopy->read(is);
        generanda_[i] = particleToCopy->copy();
        delete particleToCopy;
        generanda_[i]->setSpecies(getHandler()->getDPMBase()->speciesHandler.getObject(
                    particleToCopy_->getIndSpecies() 
                    ));
    }
    */
}

// JMFT: TODO
void PolydisperseInsertionBoundary::write(std::ostream& os) const
{
    logger(VERBOSE, "In PolydisperseInsertionBoundary::write");
    InsertionBoundary::write(os);
    os << " numberOfGeneranda " << generanda_.size() << " ";
    for (int i = 0; i < generanda_.size(); i++)
    {
        generanda_[i]->write(os);
        os << " weight " << probabilitates_[i] << " sizeDispersity " << sizeDispersities_[i];
    }
    os << " posMin " << posMin_ << " posMax " << posMax_
       << " velMin " << velMin_ << " velMax " << velMax_
       << " ";
}

std::string PolydisperseInsertionBoundary::getName() const
{
    return "PolydisperseInsertionBoundary";
}
