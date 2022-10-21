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


#include "LinearViscoelasticInteraction.h"
#include "Species/NormalForceSpecies/SPHNormalSpecies.h"
#include "Species/BaseSpecies.h"
#include "BaseInteractable.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
SPHInteraction::SPHInteraction(BaseInteractable* P, BaseInteractable* I,
                                                             unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearViscoelasticInteraction::LinearViscoelasticInteraction() finished"<<std::endl;
#endif
}

//used for mpi
SPHInteraction::SPHInteraction()
        : BaseInteraction()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearViscoelasticInteraction::LinearViscoelasticInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
SPHInteraction::SPHInteraction(const SPHInteraction& p)
        : BaseInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearViscoelasticInteraction::LinearViscoelasticInteraction(const LinearViscoelasticInteraction& p) finished"<<std::endl;
#endif
}

SPHInteraction::~SPHInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LinearViscoelasticInteraction::~LinearViscoelasticInteraction() finished"<<std::endl;
#endif
}

/*!
 * \details Calls the BaseInteraction() write function.
 * \param[in,out] os
 */
void SPHInteraction::write(std::ostream& os) const
{
    BaseInteraction::write(os);
}

/*!
 * \details Calls the BaseInteraction() read function.
 * \param[in,out] is
 */
void SPHInteraction::read(std::istream& is)
{
    BaseInteraction::read(is);
}

/*!
 * \return std::string
 */
std::string SPHInteraction::getBaseName() const
{
    return "LinearViscoelastic";
}

/*!
 * This is where the actually force law goes. In the case of SPH we are going to create an CG Point to evaluate the fields we need.
 */
void SPHInteraction::computeNormalForce()
{
    //Need add CGPoint object here
    // Compute the relative velocity vector of particle P w.r.t. I
    setRelativeVelocity(
            getP()->getVelocityAtContact(getContactPoint()) - getI()->getVelocityAtContact(getContactPoint()));
    // Compute the projection of vrel onto the normal (can be negative)
    setNormalRelativeVelocity(Vec3D::dot(getRelativeVelocity(), getNormal()));
    
    if (getOverlap() > 0) //if contact forces
    {
        const SPHNormalSpecies* species = getSpecies();
        
        Mdouble normalForce =
                species->getStiffness() * getOverlap() - species->getDissipation() * getNormalRelativeVelocity();
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
Mdouble SPHInteraction::getElasticEnergy() const
{
    if (getOverlap() > 0)
        return 0.5 * (getSpecies()->getStiffness() * mathsFunc::square(getOverlap()));
    else
        return 0.0;
}

/*!
 * \return const LinearViscoelasticNormalSpecies*
 */
const SPHNormalSpecies* SPHInteraction::getSpecies() const
{
    ///\todo you can remove this by adding the species to the Interaction instead of the base interaction
    return static_cast<const SPHNormalSpecies*>(getBaseSpecies()->getNormalForce());
}


Mdouble SPHInteraction::getElasticEnergyAtEquilibrium(Mdouble adhesiveForce) const
{
    double equilibriumOverlap = adhesiveForce / getSpecies()->getStiffness();
    return 0.5 * adhesiveForce * equilibriumOverlap;
}
