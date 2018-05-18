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

#ifndef BOUNDARIES_INSERTIONBOUNDARY_H
#define BOUNDARIES_INSERTIONBOUNDARY_H

#include "BaseBoundary.h"

class MD;
class RNG;

  /*!
   * \class InsertionBoundary
   * \brief Boundary structure for boundaries used for insertion of particles
   * \todo IFCD: Should operator= be implemented here and in the derived classes?
   */
class InsertionBoundary : public BaseBoundary
{
public:
  /*!
   * \brief Default constructor: set everything to 0/nullptr.
   */
    InsertionBoundary();
    
   /*!
    * \brief Copy constructor (with deep copy)
    */ 
    InsertionBoundary(const InsertionBoundary& other);
    
  /*!
   * \brief Destructor: delete the particle that has to be copied at every insertion.
   */
  ~InsertionBoundary() override;
    
  /*!
   * \brief Sets the particle that will be inserted and the maximum number of times for which insertion may fail.
   */
    void set(BaseParticle* particleToCopy, unsigned int maxFailed);
    
  /*!
   * \brief Virtual function that generates the intrinsic properties
   * (species, radius) of one particle. 
   * \param[in] random  Random number generator 
   */
    virtual BaseParticle* generateParticle(RNG &random)=0;

  /*!
   * \brief Purely virtual function that generates the extrinsic properties
   * (position, velocity) of a particle.
   * \param[in] p The particle to be placed
   * \param[in] random Random number generator
   * \details This should be implemented by the children such as
   * CubeInsertionBoundary, as the implementation will be geometry-dependent.
   */
    virtual void placeParticle(BaseParticle* p, RNG &random)=0;

  /*!
   * \brief Fills the boundary with particles.
   */
  void checkBoundaryBeforeTimeStep(DPMBase* md) override;

  /*!
   * \brief Gets the number of particles inserted by the boundary.
   */    
    unsigned int getNumberOfParticlesInserted() const;

    double getMassOfParticlesInserted() const;
    double getVolumeOfParticlesInserted() const;

    void reset();

    /*! 
     * \brief Turns on the InsertionBoundary.
     */
    void activate();

    /*!
     * \brief Turns off the InsertionBoundary.
     */
    void deactivate();

    
  /*!
   * \brief Sets the number of times that the wall may fail to insert a particle.
   */
    void setMaxFailed(unsigned int maxFailed);
    
  /*!
   * \brief Gets the number of times that the boundary may fail to insert a particle.
   */
    unsigned int getMaxFailed() const;
    
  /*!
   * \brief Sets the particle that will be inserted through the insertion boundary.
   */
    void setParticleToCopy(BaseParticle* particleToCopy);
    
  /*!
   * \brief Gets the particle that will be inserted through the insertion boundary.
   */
    BaseParticle* getParticleToCopy() const;
    
  /*!
   * \brief Reads the boundary's id_ and maxFailed_.
   */
    void read(std::istream& is) override;
    
  /*!
   * \brief Writes the boundary's id_ and maxFailed_.
   */
    void write(std::ostream& os) const override;

    Mdouble getVolumeFlowRate() const;

    void setVolumeFlowRate(Mdouble volumeFlowRate_);

protected:

  /*!
   * \brief Particle that will be inserted through the insertion boundary.
   */
    BaseParticle* particleToCopy_;
    
  /*!
   * \brief Number of times that the wall may fail to insert a particle.
   */
    unsigned int maxFailed_;
    
  /*!
   * \brief Number of particles that are already inserted.
   */
    unsigned int numberOfParticlesInserted_;

  /*!
   * \brief Total mass of particles inserted 
   */
    Mdouble massInserted_;

  /*!
   * \brief Total volume of particles inserted 
   */
    Mdouble volumeInserted_;

  /*! 
   * \brief The InsertionBoundary is activated by default. If the
   * InsertionBoundary is deactivated, then it introduces no particles (useful
   * for trying to maintain a certain insertion rate).
   */
    bool isActivated_;

    /*!
     * \brief defines a maximum flow rate beyond which no particles are inserted.
     */
    Mdouble volumeFlowRate_;
};

#endif
