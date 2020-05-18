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

#ifndef BOUNDARIES_POLYDISPERSEINSERTIONBOUNDARY_H
#define BOUNDARIES_POLYDISPERSEINSERTIONBOUNDARY_H

#include "BaseBoundary.h"
#include "InsertionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Math/RNG.h"
#include "Math/Vector.h"

#include <vector>

/*!
 * \class PolydisperseInsertionBoundary
 * \brief Like an InsertionBoundary but generates particles of multiple types.
 * Note that, as a child of InsertionBoundary, this class has a member called
 * particleToCopy_, which is a pointer to a particle. 
 * This pointer needs to point to something arbitrary but it doesn't matter what
 * the value is.
 */

class PolydisperseInsertionBoundary : public InsertionBoundary
{
public:
    /*!
     * \brief Constructor; sets everything to 0.
     */
    PolydisperseInsertionBoundary();
    
    /*!
     * \brief Copy constructor with deep copy.
     */
    PolydisperseInsertionBoundary(const PolydisperseInsertionBoundary& other);
    
    /*!
     * \brief Destructor: default destructor.
     */
    ~PolydisperseInsertionBoundary() override;
    
    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    PolydisperseInsertionBoundary* copy() const override;
    
    /*!
     * \brief Set position and velocity of inserted particles.
     */
    void setGeometry(int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax);
    
    /*!
     * \brief Get the particles that need to be copied 
     */
    BaseParticle* getGenerandum(unsigned int spec) const;
    
    /*!
     * \brief Add a new prototype of particle to be copied
     */
    void addGenerandum(BaseParticle* generandum, double probability, double sizeDispersity);
    
    /*!
     * \brief Change a particle to be copied
     */
    void setGenerandum(unsigned int spec, BaseParticle* generandum, double probability, double sizeDispersity);
    
    /*!
     * \brief Generates a particle from the possible species
     */
    BaseParticle* generateParticle(RNG& random) override;
    
    /*!
     * \brief Places the particle in a random position with a random velocity
     */
    void placeParticle(BaseParticle* p, RNG& random) override;
    
    /*!
     * \brief reads boundary properties from istream
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief writes boundary properties to ostream
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;

// private:
    
    /*!
     * \brief Prototypes of the particles that can be generated 
     */
    std::vector<BaseParticle*> generanda_;
    
    /*!
     * \brief The probabilities of generating each type of particle. 
     * These probabilities are not normalised.
     */
    std::vector<Mdouble> probabilitates_;
    
    /*!
     * \brief The dispersity allowed in the particle size.
     */
    std::vector<Mdouble> sizeDispersities_; // size dispersity in the radii
    
    /*!
     * \brief As in CubeInsertionBoundary. 
     * JMFT: TODO: Later we should completely separate InsertionBoundary
     * geometry from their 'intrinsics', so that the code from
     * CubeInsertionBoundary can be recycled.
     */
    Vec3D posMin_, posMax_, velMin_, velMax_;
    
    /*!
     * \brief For keeping track of how much of each prototype we have inserted.
     */
    std::vector<unsigned int> numbersInserted_;
    std::vector<Mdouble> massesInserted_;
    std::vector<Mdouble> volumesInserted_;
};

#endif
