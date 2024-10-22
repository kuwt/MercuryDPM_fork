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

#ifndef MERCURY3D_H
#define MERCURY3D_H

#include "MercuryBase.h"

class ClusterGenerator;

/*!
 * \brief This adds on the hierarchical grid code for 3D problems.
 */
class Mercury3D : public MercuryBase
{
public:

    /*!
     * \brief This is the default constructor. All it does is set sensible defaults.
     */
    Mercury3D();
    
    /*!
     * \brief Copy-constructor for creates an Mercury3D problem from an existing MD problem.
     */
    explicit Mercury3D(const DPMBase& other);
    
    /*!
     * \brief Copy-constructor.
     */
    Mercury3D(const Mercury3D& other);
    
    /*!
     * \brief Function that sets the SystemDimension and ParticleDimension to 3.
     */
    void constructor();
    
    /*!
     * \brief Returns all particles that have a contact with a given particle.
     */
    std::vector<BaseParticle*> hGridFindParticleContacts(const BaseParticle* obj) override;

protected:
    /*!
     * \brief Finds contacts between particles in the target cell.
     */
    void hGridFindContactsWithinTargetCell(int x, int y, int z, unsigned int l);
    
    /*!
     * \brief Finds contacts between the BaseParticle and the target cell.
     */
    void hGridFindContactsWithTargetCell(int x, int y, int z, unsigned int l, BaseParticle* obj);
    
    /*!
     * \brief Compute contacts with a wall.
     */
    void computeWallForces(BaseWall* w) override;


    /*!
     * \brief Finds particles within target cell and stores them in a list
     */
    void hGridFindParticlesWithTargetCell(int x, int y, int z, unsigned int l, BaseParticle* obj,
                                          std::vector<BaseParticle*>& list);
    
    /*!
     * \brief Obtains all neighbour particles of a given object, obtained from the hgrid
     */
    void hGridGetInteractingParticleList(BaseParticle* obj, std::vector<BaseParticle*>& list) override;
    
    /*!
     * \brief Finds contacts with the BaseParticle; avoids multiple checks.
     */
    void computeInternalForces(BaseParticle* obj) override;
    
    /*!
     * \brief Tests if the BaseParticle has contacts with other Particles in the target cell.
     */
    bool hGridHasContactsInTargetCell(int x, int y, int z, unsigned int l, const BaseParticle* obj) const;
    
    /*!
     * \brief Tests if a BaseParticle has any contacts in the HGrid.
     */
    bool hGridHasParticleContacts(const BaseParticle* obj) override;
    
    /*!
     * \brief Removes a BaseParticle from the HGrid.
     */
    void hGridRemoveParticle(BaseParticle* obj) override;
    
    /*!
     * \brief Updates the cell (not the level) of a BaseParticle.
     */
    void hGridUpdateParticle(BaseParticle* obj) override;

#ifdef CONTACT_LIST_HGRID
    /*!
     * \brief Adds the combination of all objects in the cell with ID (x,y,z,l) and given BaseParticle to the list of possible contacts.
     */	
    void InsertCell(int x, int y, int z, unsigned int l, BaseParticle* obj);
    
    /*!
     * \brief Add the given BaseParticle to possible interactions for all levels.
     */
    void InsertObjAgainstGrid(BaseParticle* obj);
#endif
};

#endif
