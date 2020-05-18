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


#include "HertzianSinterInteraction.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"
#include "Species/NormalForceSpecies/HertzianSinterNormalSpecies.h"
#include "Math/ExtendedMath.h"
#include <iomanip>
#include <fstream>
#include <cmath>    // std::max

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
HertzianSinterInteraction::HertzianSinterInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    maxOverlap_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianSinterInteraction::HertzianSinterInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p 
 */
HertzianSinterInteraction::HertzianSinterInteraction(const HertzianSinterInteraction& p)
        : BaseInteraction(p)
{
    maxOverlap_ = p.maxOverlap_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianSinterInteraction::HertzianSinterInteraction(const HertzianSinterInteraction &p finished"<<std::endl;
#endif
}

HertzianSinterInteraction::HertzianSinterInteraction() = default;

/*!
 *
 */
HertzianSinterInteraction::~HertzianSinterInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"HertzianSinterInteraction::~HertzianSinterInteraction() finished"<<std::endl;
#endif
}

/*!
 * \details Calls the write function of BaseInteraction().
 * \param[in,out] os
 */
void HertzianSinterInteraction::write(std::ostream& os) const
{
    BaseInteraction::write(os);
    os << " maxOverlap " << maxOverlap_;
}

/*!
 * \details Calls the read function of BaseInteraction().
 * \param[in,out] is
 */
void HertzianSinterInteraction::read(std::istream& is)
{
    BaseInteraction::read(is);
    std::string dummy;
    is >> dummy >> maxOverlap_;
}

/*!
 * \return std::string
 */
std::string HertzianSinterInteraction::getBaseName() const
{
    return "HertzianSinter";
}

/*!
 *
 */
void HertzianSinterInteraction::computeSinterForce()
{
    // Compute the relative velocity vector of particle P w.r.t. I
    setRelativeVelocity(
            getP()->getVelocityAtContact(getContactPoint()) - getI()->getVelocityAtContact(getContactPoint()));
    // Compute the projection of vrel onto the normal (can be negative)
    setNormalRelativeVelocity(Vec3D::dot(getRelativeVelocity(), getNormal()));
    
    if (getOverlap() > 0) //if contact forces
    {
        const HertzianSinterNormalSpecies* species = getSpecies();
        Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();
        
        //calculate the overlap above which the max. unloading stiffness becomes active (the 'fluid branch')
        static Mdouble maxFactor = 1
                                   - mathsFunc::square(
                cbrt((species->getLoadingModulus() + species->getCohesionModulus()) /
                     species->getUnloadingModulusMax()));
        Mdouble deltaStar = species->getPenetrationDepthMax() * effectiveDiameter / maxFactor;
        
        //increase max overlap if necessary
        if (getOverlap() > getMaxOverlap())
        {
            setMaxOverlap(getOverlap());
            std::cout << "," << getHandler()->getDPMBase()->getTime();
        }
        //limit max overlap if necessary
        if (getMaxOverlap() > deltaStar)
            setMaxOverlap(deltaStar);
        
        //calculate the unloading modulus
        Mdouble loadingCohesionModulus = species->getLoadingModulus() + species->getCohesionModulus();
        Mdouble unloadingModulus = loadingCohesionModulus
                                   + (species->getUnloadingModulusMax() - loadingCohesionModulus) *
                                     (getMaxOverlap() / deltaStar);
        
        //calculate the overlap where the force is minimal
        Mdouble factor = 1 - mathsFunc::square(cbrt(loadingCohesionModulus / unloadingModulus));
        Mdouble minOverlap = factor * maxOverlap_;
        
        //add dissipative force
        Mdouble normalForce = -species->getDissipation() * getNormalRelativeVelocity();
        
        //compute elastic force
        if (getOverlap() < minOverlap)
        {
            //decrease max overlap if in cohesive range
            std::cout << "." << getHandler()->getDPMBase()->getTime();
            setMaxOverlap(getOverlap() / factor);
        }
        else
        {
            Mdouble contactRadius = sqrt(2.0 * effectiveDiameter * (getOverlap() - minOverlap));
            normalForce += 4. / 3. * unloadingModulus * contactRadius * (getOverlap() - minOverlap);
        }
        
        setAbsoluteNormalForce(std::abs(normalForce)); //used for the friction force calculations;
        
        Mdouble contactRadius = sqrt(2.0 * effectiveDiameter * getOverlap());
        normalForce -= 4. / 3. * species->getCohesionModulus() * contactRadius * getOverlap();
        
        setForce(getNormal() * normalForce);
        setTorque(Vec3D(0.0, 0.0, 0.0));
        
        //now add the sintering model 'modified Frenkel' of the Pokula paper
        //plasticOverlap_+=species->getSinterRate()*(deltaStar-plasticOverlap_)*getHandler()->getDPMBase()->getTimeStep();
        //x/a=sqrt(2*a*del)/a
        Mdouble x = 1e-10 + sqrt(2.0 * maxOverlap_ / effectiveDiameter);
        //Mdouble x2 = x*x;
        Mdouble dx = 0.5 /
                     x;//+ x*(-0.5 + x2* (0.15625 + x2*(-0.0208333 +x2*(-0.00325521  +x2*(0.000189887 +x2*0.0000542535)))));
        Mdouble doverlap = x * dx * effectiveDiameter;
        //doverlap = 0.5/(factor*factor*plasticOverlap_);
        maxOverlap_ += species->getSinterRate() * doverlap * getHandler()->getDPMBase()->getTimeStep();
    }
    else
    {
        setAbsoluteNormalForce(0.0);
        setForce(Vec3D(0.0, 0.0, 0.0));
        setTorque(Vec3D(0.0, 0.0, 0.0));
    }
}

/*!
 *
 */
void HertzianSinterInteraction::computeNormalForce()
{
    computeSinterForce();
}

/*!
 * \return Mdouble
 */
Mdouble HertzianSinterInteraction::getElasticEnergy() const
{
    if (getOverlap() > 0)
        return 0.5 * (getSpecies()->getLoadingModulus() * mathsFunc::square(getOverlap()));
    else
        return 0.0;
    ///\todo TW this is not correct; we should count the return energy
}

/*!
 * \return const HertzianSinterNormalSpecies*
 */
const HertzianSinterNormalSpecies* HertzianSinterInteraction::getSpecies() const
{
    return static_cast<const HertzianSinterNormalSpecies*>(getBaseSpecies()->getNormalForce());
}

/*!
 * \return Mdouble plasticOverlap_
 */
Mdouble HertzianSinterInteraction::getMaxOverlap() const
{
    return maxOverlap_;
}

/*!
 * \param[in] maxOverlap
 */
void HertzianSinterInteraction::setMaxOverlap(const Mdouble maxOverlap)
{
    maxOverlap_ = maxOverlap;
}

/*!
 * \return Mdouble
 */
Mdouble HertzianSinterInteraction::getUnloadingModulus() const
{
    const HertzianSinterNormalSpecies* species = getSpecies();
    Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();
    Mdouble deltaMaxFluid = species->getPenetrationDepthMax() * effectiveDiameter /
                            (1.0 - species->getLoadingModulus() / species->getUnloadingModulusMax());
    if (getOverlap() > deltaMaxFluid)
        return species->getUnloadingModulusMax();
    else
        return species->getLoadingModulus() +
               (species->getUnloadingModulusMax() - species->getLoadingModulus()) * getMaxOverlap() / deltaMaxFluid;
}

