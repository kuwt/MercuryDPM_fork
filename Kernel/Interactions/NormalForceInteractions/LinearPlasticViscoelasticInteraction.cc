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


#include "LinearPlasticViscoelasticInteraction.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"
#include <iomanip>
#include <fstream>
#include <Species/NormalForceSpecies/LinearPlasticViscoelasticNormalSpecies.h>
#include <cmath>    // std::max

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
LinearPlasticViscoelasticInteraction::LinearPlasticViscoelasticInteraction(BaseInteractable* P, BaseInteractable* I,
                                                                           unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    maxOverlap_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteraction::LinearPlasticViscoelasticInteraction() finished"<<std::endl;
#endif
}

//used for mpi
LinearPlasticViscoelasticInteraction::LinearPlasticViscoelasticInteraction()
        : BaseInteraction()
{
    maxOverlap_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteraction::LinearPlasticViscoelasticInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p 
 */
LinearPlasticViscoelasticInteraction::LinearPlasticViscoelasticInteraction(
        const LinearPlasticViscoelasticInteraction& p)
        : BaseInteraction(p)
{
    maxOverlap_ = p.maxOverlap_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteraction::LinearPlasticViscoelasticInteraction(const LinearPlasticViscoelasticInteraction &p finished"<<std::endl;
#endif
}

/*!
 *
 */
LinearPlasticViscoelasticInteraction::~LinearPlasticViscoelasticInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteraction::~LinearPlasticViscoelasticInteraction() finished"<<std::endl;
#endif
}

/*!
 * \details Calls the write function of BaseInteraction().
 * \param[in,out] os
 */
void LinearPlasticViscoelasticInteraction::write(std::ostream& os) const
{
    BaseInteraction::write(os);
    os << " maxOverlap " << maxOverlap_;
}

/*!
 * \details Calls the read function of BaseInteraction().
 * \param[in,out] is
 */
void LinearPlasticViscoelasticInteraction::read(std::istream& is)
{
    BaseInteraction::read(is);
    std::string dummy;
    is >> dummy >> maxOverlap_;
    //helpers::readOptionalVariable<Mdouble>(is, "maxOverlap", maxOverlap_);
}

/*!
 * \return std::string
 */
std::string LinearPlasticViscoelasticInteraction::getBaseName() const
{
    return "LinearPlasticViscoelastic";
}

/*!
 *
 */
void LinearPlasticViscoelasticInteraction::computeNormalForce()
{
    // Compute the relative velocity vector of particle P w.r.t. I
    setRelativeVelocity(
            getP()->getVelocityAtContact(getContactPoint()) - getI()->getVelocityAtContact(getContactPoint()));
    // Compute the projection of vrel onto the normal (can be negative)
    setNormalRelativeVelocity(Vec3D::dot(getRelativeVelocity(), getNormal()));

    if (getOverlap() > 0) //if contact forces
    {
        const LinearPlasticViscoelasticNormalSpecies* species = getSpecies();

        //calculate the overlap above which the max. unloading stiffness becomes active (the 'fluid branch')
        // Modified by Paolo
        // Previously:
        // const Mdouble effectiveDiameter = species->getConstantRestitution()?1.0:(2.0 * getEffectiveRadius());
        // This has been modified so that, independently on whether or not the users sets constant restitution,
        // this parameter is set to its relative value, i.e. the penetration depth max divided by the radius of a single particle.
        const Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();
        const Mdouble deltaStar = (species->getUnloadingStiffnessMax()
                                   / (species->getUnloadingStiffnessMax() - species->getLoadingStiffness()))
                                  * species->getPenetrationDepthMax() * effectiveDiameter;

        //increase max overlap if necessary
        if (getOverlap() > getMaxOverlap())
            setMaxOverlap(std::min(deltaStar, getOverlap()));

        //calculate the unloading stiffness
        const Mdouble unloadingStiffness =
                species->getDoConstantUnloadingStiffness()?species->getUnloadingStiffnessMax():
                (species->getLoadingStiffness() +
                 (species->getUnloadingStiffnessMax() - species->getLoadingStiffness()) * (getMaxOverlap() / deltaStar));

        //calculate the overlap where the force is zero
        const Mdouble equilibriumOverlap =
                (unloadingStiffness - species->getLoadingStiffness()) / unloadingStiffness * maxOverlap_;

        //compute elastic force
        Mdouble normalForce = unloadingStiffness * (getOverlap() - equilibriumOverlap);

        //decrease max overlap if necessary
        const Mdouble cohesiveForce = -species->getCohesionStiffness() * getOverlap();
        if (normalForce < cohesiveForce)
        {
            setMaxOverlap((unloadingStiffness + species->getCohesionStiffness())
                          / (unloadingStiffness - species->getLoadingStiffness()) * getOverlap());
            //only necessary because the timeStep is finite:
            normalForce = cohesiveForce;
        }

        //add dissipative force
        normalForce -= species->getDissipation() * getNormalRelativeVelocity();

        if (species->getConstantRestitution()) normalForce *= 2.0*getEffectiveMass();
        setAbsoluteNormalForce(std::abs(normalForce)); //used for further corce calculations;
        setForce(getNormal() * normalForce);
        setTorque(Vec3D(0.0, 0.0, 0.0));
    }
    else
    {
        setAbsoluteNormalForce(0.0);
        setForce(Vec3D(0.0, 0.0, 0.0));
        setTorque(Vec3D(0.0, 0.0, 0.0));
    }
}

/*!
 * \return Mdouble
 */
Mdouble LinearPlasticViscoelasticInteraction::getElasticEnergy() const
{
    Mdouble energy = getOverlap() > 0 ? 0.5 * (getSpecies()->getLoadingStiffness() * mathsFunc::square(getOverlap())) : 0.0;
    if (getBaseSpecies()->getNormalForce()->getConstantRestitution()) energy *= 2.0*getEffectiveMass();
    return energy;
    ///\todo TW this is not correct; we should count the return energy
}

/*!
 * \return const LinearPlasticViscoElasticNormalSpecies*
 */
const LinearPlasticViscoelasticNormalSpecies* LinearPlasticViscoelasticInteraction::getSpecies() const
{
    return static_cast<const LinearPlasticViscoelasticNormalSpecies*>(getBaseSpecies()->getNormalForce());
;
}

/*!
 * \return Mdouble plasticOverlap_
 */
Mdouble LinearPlasticViscoelasticInteraction::getMaxOverlap() const
{
    return maxOverlap_;
}

/*!
 * \param[in] maxOverlap
 */
void LinearPlasticViscoelasticInteraction::setMaxOverlap(const Mdouble maxOverlap)
{
    maxOverlap_ = maxOverlap;
}

/*!
 * \return Mdouble
 */
Mdouble LinearPlasticViscoelasticInteraction::getUnloadingStiffness() const
{
    const LinearPlasticViscoelasticNormalSpecies* species = getSpecies();
    // Modified by Paolo: read comment in computeNormalForce().
    // previously: const Mdouble effectiveDiameter = species->getConstantRestitution()?1.0:(2.0 * getEffectiveRadius());
    const Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();
    Mdouble deltaMaxFluid = species->getPenetrationDepthMax() * effectiveDiameter /
                            (1.0 - species->getLoadingStiffness() / species->getUnloadingStiffnessMax());
    if (getOverlap() > deltaMaxFluid)
        return species->getUnloadingStiffnessMax();
    else
        return species->getLoadingStiffness() +
               (species->getUnloadingStiffnessMax() - species->getLoadingStiffness()) * getMaxOverlap() / deltaMaxFluid;
}

