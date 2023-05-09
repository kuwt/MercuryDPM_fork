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


#include "LinearPlasticViscoelasticInteractionCutoff.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"
#include <iomanip>
#include <fstream>
#include "LinearPlasticViscoelasticNormalSpeciesExtended.h"
#include <cmath>    // std::max
/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
LinearPlasticViscoelasticInteractionCutoff::LinearPlasticViscoelasticInteractionCutoff(BaseInteractable* P, BaseInteractable* I, Mdouble timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    maxOverlap_=0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteractionCutoff::LinearPlasticViscoelasticInteractionCutoff() finished"<<std::endl;
#endif
}
/*!
 * \param[in] p
 */
LinearPlasticViscoelasticInteractionCutoff::LinearPlasticViscoelasticInteractionCutoff(const LinearPlasticViscoelasticInteractionCutoff &p)
        : BaseInteraction(p)
{
    maxOverlap_=p.maxOverlap_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteractionCutoff::LinearPlasticViscoelasticInteractionCutoff(const LinearPlasticViscoelasticInteractionCutoff &p finished"<<std::endl;
#endif
}
/*!
 *
 */
LinearPlasticViscoelasticInteractionCutoff::~LinearPlasticViscoelasticInteractionCutoff()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteractionCutoff::~LinearPlasticViscoelasticInteractionCutoff() finished"<<std::endl;
#endif
}

/*!
 * \details Calls the write function of BaseInteraction().
 * \param[in,out] os
 */
void LinearPlasticViscoelasticInteractionCutoff::write(std::ostream& os) const
{
    BaseInteraction::write(os);
    os << " maxOverlap " << maxOverlap_;
}
/*!
 * \details Calls the read function of BaseInteraction().
 * \param[in,out] is
 */
void LinearPlasticViscoelasticInteractionCutoff::read(std::istream& is)
{
    BaseInteraction::read(is);
    helpers::readOptionalVariable<Mdouble>(is,"maxOverlap",maxOverlap_);
}
/*!
 * \return std::string
 */
std::string LinearPlasticViscoelasticInteractionCutoff::getBaseName() const
{
    return "LinearPlasticViscoelastic";
}
/*!
 *
 */
void LinearPlasticViscoelasticInteractionCutoff::computeLinearPlasticViscoelasticForce()
{
    // Compute the relative velocity vector of particle P w.r.t. I
    setRelativeVelocity(getP()->getVelocityAtContact(getContactPoint()) - getI()->getVelocityAtContact(getContactPoint()));
    // Compute the projection of vrel onto the normal (can be negative)
    setNormalRelativeVelocity(Vec3D::dot(getRelativeVelocity(), getNormal()));

    if (getOverlap() > 0) //if contact forces
    {
        const LinearPlasticViscoelasticNormalSpeciesExtended* species = getSpecies();

        //calculate the overlap above which the max. unloading stiffness becomes active (the 'fluid branch')
        Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();
        Mdouble deltaStar = (species->getUnloadingStiffnessMax()
            / (species->getUnloadingStiffnessMax() - species->getLoadingStiffness()))
            * species->getPenetrationDepthMax() * effectiveDiameter;

        //increase max overlap if necessary
        if (getOverlap()>getMaxOverlap())
            setMaxOverlap(std::min(deltaStar, getOverlap()));

        //calculate the unloading stiffness
        Mdouble unloadingStiffness = species->getLoadingStiffness()
            + (species->getUnloadingStiffnessMax() - species->getLoadingStiffness()) * (getMaxOverlap() / deltaStar);

        //calculate the overlap where the force is zero
        Mdouble equilibriumOverlap = (unloadingStiffness - species->getLoadingStiffness()) / unloadingStiffness * maxOverlap_;

        //compute elastic force
        Mdouble normalForce = unloadingStiffness * (getOverlap() - equilibriumOverlap);

        //decrease max overlap if necessary
        Mdouble cohesiveForce = -species->getCohesionStiffness() * getOverlap();
        if (normalForce < cohesiveForce) {
            setMaxOverlap((unloadingStiffness + species->getCohesionStiffness())
                / (unloadingStiffness - species->getLoadingStiffness()) * getOverlap());
            //only necessary because the timeStep is finite:
            normalForce = cohesiveForce;
        }

        //calculate delta min
        Mdouble deltaMin = (unloadingStiffness - species->getLoadingStiffness()) / (unloadingStiffness + species->getCohesionStiffness()) * maxOverlap_;

        //resets the history by setting the history parameter maxOverlap_ to 0
        // if (getOverlap() < deltaMin*species->getEpsilon())
        if (normalForce < 0.0 && getOverlap() < species->getEpsilon()*deltaMin)
        {
           // std::cout << "epsilon = " << species->getEpsilon() << ", dMin = " << deltaMin << ", dEq = " << equilibriumOverlap << ", d = " << getOverlap() << ", oldF = " << normalForce;
           normalForce = species->getLoadingStiffness() * getOverlap();
           // setMaxOverlap(getOverlap());   // NEVER PUT THIS BACK IN!
           // std::cout << ", newF = " << normalForce << std::endl;
        }


        //add dissipative force
        normalForce -= species->getDissipation() * getNormalRelativeVelocity();

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
 *
 */
void LinearPlasticViscoelasticInteractionCutoff::computeNormalForce()
{
    computeLinearPlasticViscoelasticForce();
}
/*!
 * \return Mdouble
 */
Mdouble LinearPlasticViscoelasticInteractionCutoff::getElasticEnergy() const
{
   if (getOverlap() > 0)
        return 0.5 * (getSpecies()->getLoadingStiffness() * mathsFunc::square(getOverlap()));
    else
        return 0.0;
  ///\todo TW this is not correct; we should count the return energy
}
/*!
 * \return const LinearPlasticViscoElasticNormalSpecies*
 */
const LinearPlasticViscoelasticNormalSpeciesExtended* LinearPlasticViscoelasticInteractionCutoff::getSpecies() const
{
    return dynamic_cast<const LinearPlasticViscoelasticNormalSpeciesExtended*>(getBaseSpecies());
}
/*!
 * \return Mdouble plasticOverlap_
 */
Mdouble LinearPlasticViscoelasticInteractionCutoff::getMaxOverlap() const
{
    return maxOverlap_;
}
/*!
 * \param[in] maxOverlap
 */
void LinearPlasticViscoelasticInteractionCutoff::setMaxOverlap(const Mdouble maxOverlap)
{
    maxOverlap_ = maxOverlap;
}


/*!
 * \return Mdouble
 */
Mdouble LinearPlasticViscoelasticInteractionCutoff::getUnloadingStiffness() const
{
    const LinearPlasticViscoelasticNormalSpeciesExtended* species = getSpecies();
    Mdouble effectiveDiameter = 2.0*getEffectiveRadius();
    Mdouble deltaMaxFluid = species->getPenetrationDepthMax() * effectiveDiameter / (1.0-species->getLoadingStiffness()/species->getUnloadingStiffnessMax());
    if (getOverlap() > deltaMaxFluid)
        return species->getUnloadingStiffnessMax();
    else
        return species->getLoadingStiffness() + (species->getUnloadingStiffnessMax() - species->getLoadingStiffness()) * getMaxOverlap()/deltaMaxFluid;
}
