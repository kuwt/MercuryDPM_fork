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

#include <random>
#include "BaseClusterInsertionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Math/RNG.h"

/*!
 * \details Default constructor;
 */
BaseClusterInsertionBoundary::BaseClusterInsertionBoundary() : InsertionBoundary()
{
    /*
     * CLUSTER VALUES
     */

    clusterSpecies_ = new LinearPlasticViscoelasticFrictionSpecies;

    //Position
    position_ = {0,0,0};

    nClusterInserted_ = 0;

    // Time
    collisionTimeOverTimeStep_ = 50;

    // Particles
    sizeDispersityParticle_ = 1;
    radiusParticle_ = 0;
    nParticles_ = 0;


    // Central force
    velocityDampingModulus_ = 0.9;

    // Data analysis
    nInternalStructurePoints_ = 1000000;

    // Energy
    energyRatioTolerance_ = 1.0e-7;

    // File
    isCdatOutputOn_ = false;
    isOverlOutputOn_ = false;
    isAmatOutputOn_ = false;
    isIntStrucOutputOn_ = false;
    isVtkOutputOn_ = false;
    isRestartOutputOn_ = false;
    isFStatOutputOn_ = false;
    isEneOutputOn_ = false;

    radiusParticle_ = 0;


    /*
     * TYPICAL INSERTION BOUNDARY VALUES
     */

    posMin_ = Vec3D(0.0, 0.0, 0.0);
    posMax_ = Vec3D(0.0, 0.0, 0.0);
    velMin_ = Vec3D(0.0, 0.0, 0.0);
    velMax_ = Vec3D(0.0, 0.0, 0.0);

    setRadiusParticleAndNotNumberOfParticles_ = true;

    randomised_ = true;

    logger(DEBUG, "BaseClusterInsertionBoundary::BaseClusterInsertionBoundary() finished");

}

/*!
 * \details Copy constructor
 */
BaseClusterInsertionBoundary::BaseClusterInsertionBoundary(const BaseClusterInsertionBoundary& other)
        : InsertionBoundary(other)
{
    /*
     * CLUSTER VALUES
     */

    clusterSpecies_ = other.clusterSpecies_;

    //Position
    position_ = other.position_;

    nClusterInserted_ = other.nClusterInserted_;

    // Time
    collisionTimeOverTimeStep_ = other.collisionTimeOverTimeStep_;

    // Particles
    sizeDispersityParticle_ = other.sizeDispersityParticle_;
    radiusParticle_ = other.radiusParticle_;
    nParticles_ = other.nParticles_;

    nParticles_ = other.nParticles_;


    // Central force
    velocityDampingModulus_ = other.velocityDampingModulus_;

    // Data analysis
    nInternalStructurePoints_ = other.nInternalStructurePoints_;

    // Energy
    energyRatioTolerance_ = other.energyRatioTolerance_;

    // File
    isCdatOutputOn_ = other.isCdatOutputOn_;
    isOverlOutputOn_ = other.isOverlOutputOn_;
    isAmatOutputOn_ = other.isAmatOutputOn_;
    isIntStrucOutputOn_ = other.isIntStrucOutputOn_;
    isVtkOutputOn_ = other.isVtkOutputOn_;
    isRestartOutputOn_ = other.isRestartOutputOn_;
    isFStatOutputOn_ = other.isFStatOutputOn_;
    isEneOutputOn_ = other.isEneOutputOn_;



    radiusParticle_ = other.radiusParticle_;
    /*
     * TYPICAL INSERTION BOUNDARY VALUES
     */

    posMin_ = other.posMin_;
    posMax_ = other.posMax_;
    velMin_ = other.velMin_;
    velMax_ = other.velMax_;

    setRadiusParticleAndNotNumberOfParticles_ = other.setRadiusParticleAndNotNumberOfParticles_;
    
    randomised_=other.randomised_;

}

/*!
 * \details Destructor. Since there are no pointers in this class, there is no
 *          need for any actions here.
 */
BaseClusterInsertionBoundary::~BaseClusterInsertionBoundary()
= default;

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer.
 * \return      pointer to the copy on the heap
 */
BaseClusterInsertionBoundary* BaseClusterInsertionBoundary::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "BaseClusterInsertionBoundary::copy() const finished" << std::endl;
#endif
    return new BaseClusterInsertionBoundary(*this);
}

//!\brief this turns off the randomise(): created for UnitTests.
void BaseClusterInsertionBoundary::setRandomised(bool randomised){
    randomised_ = randomised;
}

bool BaseClusterInsertionBoundary::getRandomised(){
    return randomised_;
}

unsigned int BaseClusterInsertionBoundary::getNumberOfClusterInserted(){
    return nClusterInserted_;
}

void BaseClusterInsertionBoundary::setRadiusMicroParticle(Mdouble radiusMicroParticle)
{
    if (radiusMicroParticle <= 0)
        logger(ERROR, "radiusParticle must be greater than zero. radiusMicroParticle = %", radiusMicroParticle);
    else {
        radiusParticle_ = radiusMicroParticle;
        setRadiusParticleAndNotNumberOfParticles_ = true;
    }
}

void BaseClusterInsertionBoundary::setRadiusRange(Mdouble radMin, Mdouble radMax)
{
    radMin_ = radMin;
    radMax_ = radMax;
}

void BaseClusterInsertionBoundary::setGeometry(Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax)
{
    posMin_ = posMin;
    posMax_ = posMax;
    velMin_ = velMin;
    velMax_ = velMax;
}

void BaseClusterInsertionBoundary::setVelocityRange(Vec3D velMin, Vec3D velMax)
{
    velMin_ = velMin;
    velMax_ = velMax;
}

void BaseClusterInsertionBoundary::setAdditionalClusterProperties(Mdouble collisionTimeOverTimeStep, Mdouble velocityDampingModulus, Mdouble energyRatioTolerance)
{
    collisionTimeOverTimeStep_ = collisionTimeOverTimeStep;
    velocityDampingModulus_ = velocityDampingModulus;
    energyRatioTolerance_ = energyRatioTolerance;

}

void BaseClusterInsertionBoundary::setOutputClusterProperties(bool doCdatOutput, bool doOverlOutput, bool doAmatOutput, bool doIntStrucOutput,
                                                              bool doVtkOutput, bool doRestartOutput, bool doFStatOutput, bool doEneOutput)
{
    isCdatOutputOn_ = doCdatOutput;
    isOverlOutputOn_ = doOverlOutput;
    isAmatOutputOn_ = doAmatOutput;
    isIntStrucOutputOn_ = doIntStrucOutput;
    isVtkOutputOn_ = doVtkOutput;
    isRestartOutputOn_ = doRestartOutput;
    isFStatOutputOn_ = doFStatOutput;
    isEneOutputOn_ = doEneOutput;

}

//!\brief function overridden in children classes FixedClusterInsertionBoundary and RandomClusterInsertionBoundary
void BaseClusterInsertionBoundary::placeParticle(BaseParticle* p, RNG& random)
{
}

/*!
 * \details Is used to fill the insides of the boundary with clusters until
 * it is filled up. Function overridden in children classes FixedClusterInsertionBoundary
 *  and RandomClusterInsertionBoundary
 * \param[in,out] md    the problem's DPMBase object
 * \todo rename to something like "insertUntilMaxFailed"?
 */
void BaseClusterInsertionBoundary::checkBoundaryBeforeTimeStep(DPMBase* md)
{
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in,out] is        the istream
 */
void BaseClusterInsertionBoundary::read(std::istream& is)
{
    InsertionBoundary::read(is);
    std::string dummy;

    is >> dummy >> posMin_
       >> dummy >> posMax_;
    is >> dummy >> velMin_
       >> dummy >> velMax_;
    is >> dummy >> radMin_;
    is >> dummy >> radMax_;
    
    is >> dummy >> nClusterInserted_
         >> dummy >> radiusParticle_;
    is >> dummy >> sizeDispersityParticle_
         >> dummy >> velocityDampingModulus_
         >> dummy >> nInternalStructurePoints_;
    is >> dummy >> energyRatioTolerance_
         >> dummy >> collisionTimeOverTimeStep_
         >> dummy >> setRadiusParticleAndNotNumberOfParticles_;
    is >> dummy >> isCdatOutputOn_
         >> dummy >> isOverlOutputOn_
         >> dummy >> isAmatOutputOn_
         >> dummy >> isIntStrucOutputOn_
         >> dummy >> isVtkOutputOn_
         >> dummy >> isRestartOutputOn_
         >> dummy >> isFStatOutputOn_
         >> dummy >> isEneOutputOn_;

}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream
 */
void BaseClusterInsertionBoundary::write(std::ostream& os) const
{
    InsertionBoundary::write(os);
    os << " posMin " << posMin_
       << " posMax " << posMax_
       << " velMin " << velMin_
       << " velMax " << velMax_;
    os << " radMin " << radMin_
       << " radMax " << radMax_;

    os << " nClusterInserted " << nClusterInserted_ <<
          " radiusParticle " << radiusParticle_ <<
          " sizeDispersityParticle " << sizeDispersityParticle_ <<
          " velocityDampingModulus " << velocityDampingModulus_ <<
          " nInternalStructurePoints " << nInternalStructurePoints_ <<
          " energyRatioTolerance " << energyRatioTolerance_ <<
          " collisionTimeOverTimeStep " << collisionTimeOverTimeStep_ <<
          " setRadiusParticle " << setRadiusParticleAndNotNumberOfParticles_ <<
          " isCdatOutputOn " << isCdatOutputOn_ <<
          " isOverlOutputOn " << isOverlOutputOn_ <<
          " isAmatOutputOn " << isAmatOutputOn_ <<
          " isIntStrucOutputOn " << isIntStrucOutputOn_ <<
          " isVtkOutputOn " << isVtkOutputOn_ <<
          " isRestartOutputOn " << isRestartOutputOn_ <<
          " isFStatOutputOn " << isFStatOutputOn_ <<
          " isEneOutputOn " << isEneOutputOn_;

}








/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'ClusterInsertionBoundary'
 */
std::string BaseClusterInsertionBoundary::getName() const
{
    return "BaseClusterInsertionBoundary";
}