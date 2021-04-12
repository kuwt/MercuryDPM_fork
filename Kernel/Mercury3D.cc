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

#include "Mercury3D.h"
#include "Particles/BaseParticle.h"

Mercury3D::Mercury3D()
{
    constructor();
    logger(DEBUG, "Mercury3D::Mercury3D() finished");
}

/*!
 * \param[in] other Mercury3D that must be copied.
 * \details Copy constructor, note that the copy-constructor of DPMBase has to 
 * be called because the link from DPMBase to MercuryBase is virtual.
 */
Mercury3D::Mercury3D(const Mercury3D& other)
        : DPMBase(other), MercuryBase(other)
{
    logger(DEBUG, "Mercury3D::Mercury3D(Mercury3D& other) copy constructor finished.");
}

/*!
 * \param[in] other DPMBase which has to be copied and converted to a Mercury3D.
 * \details Constructor that makes a Mercury3D out of a DPMBase.
 * The "copy"-constructor of DPMBase has to be called because the link 
 * from DPMBase to MercuryBase is virtual.
 */
Mercury3D::Mercury3D(const DPMBase& other)
        : DPMBase(other), MercuryBase()
{
    constructor();
    logger(DEBUG, "Mercury3D::Mercury3D(DPMBase& other) constructor finished");
}

void Mercury3D::constructor()
{
    setParticleDimensions(3);
    setSystemDimensions(3);
}

/*!
 * \param[in] x Coordinate of the target cell in x direction.
 * \param[in] y Coordinate of the target cell in y direction.
 * \param[in] z Coordinate of the target cell in z direction.
 * \param[in] l Level in the HGrid of the target cell.
 * \details Computes all collision between particles in the same bucket as cell 
 * (x,y,z,l), please note that all the particles are in the same cell.
 */
void Mercury3D::hGridFindContactsWithinTargetCell(int x, int y, int z, unsigned int l)
{
    HGrid* const hgrid = getHGrid();
    const unsigned int bucket = hgrid->computeHashBucketIndex(x, y, z, l);
    
    ///\todo replace this generic check of the each bucket to checking only the object to avoid the critical
    //Check if this function is already applied to this bucket
    bool bucketIsChecked;
    #pragma omp critical
    {
        bucketIsChecked = hgrid->getBucketIsChecked(bucket);
        hgrid->setBucketIsChecked(bucket);
    }
    if (bucketIsChecked) return;
    
    BaseParticle* p1 = hgrid->getFirstBaseParticleInBucket(bucket);
    while (p1 != nullptr)
    {
        BaseParticle* p2 = p1->getHGridNextObject();
        while (p2 != nullptr)
        {
            ///\bug TW: This check is not necessary, I believe. This is the most-expensive function in most codes (the two checks in this function slows down granular jet by 15%) and the selftests are not affected.
            ///\bug DK: I do think this is necessary, for example: If two cells hash to the same bucket and a particle in one of these cells check for collisions with the other cell. Then due to the hashing collision it also gets all particles in it's own cell and thus generating false collisions.
            //Check if the BaseParticle* p1 and BaseParticle* p2 are really in the same cell (i.e. no hashing error has occurred)
            if (p1->getHGridCell() == (p2->getHGridCell()))
            {
                computeInternalForce(p1, p2);
            }
            p2 = p2->getHGridNextObject();
        }
        p1 = p1->getHGridNextObject();
    }
}

/*!
 * \param[in] x     The coordinate of the target cell in x direction.
 * \param[in] y     The coordinate of the target cell in y direction.
 * \param[in] z     The coordinate of the target cell in z direction.
 * \param[in] l     The level in the HGrid of the target cell.
 * \param[in] obj   A pointer to the BaseParticle for which we want to have interactions.
 * \details Computes all collisions between given BaseParticle and particles in 
 * cell (x,y,z,l). This is done by first checking if the BaseParticle is indeed from
 * another cell, then for all BaseParticle in the target cell it is checked what
 * the forces between that BaseParticle and given BaseParticle are.
 */
void Mercury3D::hGridFindContactsWithTargetCell(int x, int y, int z, unsigned int l, BaseParticle* const obj)
{
    //Check if the object is not in the same cell as being checked, CheckCell_current should handle these cases.
    //TW a speedcheck revealed that this check costs a 10% performance decrease; it's only a safety check, so I made it an assert.
    logger.assert(!obj->getHGridCell().equals(x, y, z, l),
                  "hGridFindContactsWithTargetCell should not be called if object is in the same cell");
    
    HGrid* const hgrid = getHGrid();
    
    // Calculate the bucket
    const unsigned int bucket = hgrid->computeHashBucketIndex(x, y, z, l);
    
    // Loop through all objects in the bucket to find nearby objects
    for (BaseParticle* p = hgrid->getFirstBaseParticleInBucket(bucket); p != nullptr; p = p->getHGridNextObject())
    {
        //This is the most-expensive function in most codes (the two checks in this function slows down granular jet by 15%). It is neccesary, for example: If two cells hash to the same bucket and a particle in one of these cells check for collisions with the other cell. Then due to the hashing collision it also gets all particles in it's own cell and thus generating false collisions.
        //Check if the BaseParticle *p really is in the target cell (i.e. no hashing error has occurred)
        //TW speedcheck revealed that this pre-check is cheaper than allowing computeInternalForces to sort out mismatches; even if a large number of hash cells (10*Np) is used.
        if (p->getHGridCell().equals(x, y, z, l))
        {
            if (Vec3D::getDistanceSquared(p->getPosition(),obj->getPosition()) < mathsFunc::square(p->getMaxInteractionRadius()+obj->getMaxInteractionRadius()))
                computeInternalForce(obj, p);
        }
    }
}

//TODO document
//TODO obj kan weg?
void Mercury3D::hGridFindParticlesWithTargetCell(int x, int y, int z, unsigned int l, BaseParticle* obj,
                                                 std::vector<BaseParticle*>& list)
{
    HGrid* const hgrid = getHGrid();
    
    // Calculate the bucket
    const unsigned int bucket = hgrid->computeHashBucketIndex(x, y, z, l);
    
    // Loop through all objects in the bucket to find nearby objects
    BaseParticle* p = hgrid->getFirstBaseParticleInBucket(bucket);
    while (p != nullptr)
    {
        if (p->getHGridCell().equals(x, y, z, l))
        {
            list.push_back(p);
        }
        p = p->getHGridNextObject();
    }
}

void Mercury3D::hGridGetInteractingParticleList(BaseParticle* obj, std::vector<BaseParticle*>& list)
{
    HGrid* hgrid = getHGrid();
    
    ///\bug find out why this is necessary; if this is not there, the code sometimes segfaults.
    if (hGridNeedsRebuilding())
    {
        hGridRebuild();
        hgrid = getHGrid();
    }
    logger(DEBUG, "hgrid %, object %", hgrid, obj);
    int occupiedLevelsMask = hgrid->getOccupiedLevelsMask() >> obj->getHGridLevel();
    for (unsigned int level = 0; level < hgrid->getNumberOfLevels(); level++)
    {
        // If no objects in rest of grid, stop now
        if (occupiedLevelsMask == 0)
        {
            break;
        }
        
        // If no objects at this level, go on to the next level
        if ((occupiedLevelsMask & 1) == 0)
        {
            continue;
        }
        
        const Mdouble inv_size = hgrid->getInvCellSize(level);
        const int xs = static_cast<int>(std::floor(
                (obj->getPosition().X - obj->getMaxInteractionRadius()) * inv_size - 0.5));
        const int xe = static_cast<int>(std::floor(
                (obj->getPosition().X + obj->getMaxInteractionRadius()) * inv_size + 0.5));
        const int ys = static_cast<int>(std::floor(
                (obj->getPosition().Y - obj->getMaxInteractionRadius()) * inv_size - 0.5));
        const int ye = static_cast<int>(std::floor(
                (obj->getPosition().Y + obj->getMaxInteractionRadius()) * inv_size + 0.5));
        const int zs = static_cast<int>(std::floor(
                (obj->getPosition().Z - obj->getMaxInteractionRadius()) * inv_size - 0.5));
        const int ze = static_cast<int>(std::floor(
                (obj->getPosition().Z + obj->getMaxInteractionRadius()) * inv_size + 0.5));
        for (int x = xs; x <= xe; ++x)
        {
            for (int y = ys; y <= ye; ++y)
            {
                for (int z = zs; z <= ze; ++z)
                {
                    hGridFindParticlesWithTargetCell(x, y, z, level, obj, list);
                }
            }
        }
    }
}

/*!
 * \param[in] obj   A pointer to the BaseParticle for which we want to check for contacts.
 * \details         Computes all collision between given BaseParticle and all other 
 *                  particles in the grid. Please note that we're looking only one way, so that 
 *                  interactions are not detected twice.
 */
void Mercury3D::computeInternalForces(BaseParticle* obj)
{
    HGrid* const hgrid = getHGrid();
    const unsigned int startLevel = obj->getHGridLevel();
    
    if (getHGridMethod() == TOPDOWN)
    {
        int occupiedLevelsMask = hgrid->getOccupiedLevelsMask();
        for (unsigned int level = 0; level <= startLevel && occupiedLevelsMask != 0; occupiedLevelsMask >>= 1, level++)
        {
            // If no objects at this level, go on to the next level
            if ((occupiedLevelsMask & 1) == 0)
            {
                continue;
            }
            
            if (level == startLevel)
            {
                const int x = obj->getHGridX();
                const int y = obj->getHGridY();
                const int z = obj->getHGridZ();
                
                hGridFindContactsWithinTargetCell(x, y, z, level);
                hGridFindContactsWithTargetCell(x + 1, y - 1, z, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y, z, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y + 1, z, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y - 1, z + 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y, z + 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y + 1, z + 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y - 1, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y + 1, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x, y + 1, z, level, obj);
                hGridFindContactsWithTargetCell(x, y, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x, y + 1, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x, y + 1, z + 1, level, obj);
            }
            else
            {
                const Mdouble inv_size = getHGrid()->getInvCellSize(level);
                const int xs = static_cast<int>(std::floor(
                        (obj->getPosition().X - obj->getMaxInteractionRadius()) * inv_size - 0.5));
                const int xe = static_cast<int>(std::floor(
                        (obj->getPosition().X + obj->getMaxInteractionRadius()) * inv_size + 0.5));
                const int ys = static_cast<int>(std::floor(
                        (obj->getPosition().Y - obj->getMaxInteractionRadius()) * inv_size - 0.5));
                const int ye = static_cast<int>(std::floor(
                        (obj->getPosition().Y + obj->getMaxInteractionRadius()) * inv_size + 0.5));
                const int zs = static_cast<int>(std::floor(
                        (obj->getPosition().Z - obj->getMaxInteractionRadius()) * inv_size - 0.5));
                const int ze = static_cast<int>(std::floor(
                        (obj->getPosition().Z + obj->getMaxInteractionRadius()) * inv_size + 0.5));
                for (int x = xs; x <= xe; ++x)
                {
                    for (int y = ys; y <= ye; ++y)
                    {
                        for (int z = zs; z <= ze; ++z)
                        {
                            hGridFindContactsWithTargetCell(x, y, z, level, obj);
                        }
                    }
                }
            }
        }
    }
    else
    {
        int occupiedLevelsMask = hgrid->getOccupiedLevelsMask() >> obj->getHGridLevel();
        for (unsigned int level = startLevel; level < hgrid->getNumberOfLevels(); occupiedLevelsMask >>= 1, level++)
        {
            // If no objects in rest of grid, stop now
            if (occupiedLevelsMask == 0)
            {
                break;
            }
            
            // If no objects at this level, go on to the next level
            if ((occupiedLevelsMask & 1) == 0)
            {
                continue;
            }
            
            if (level == startLevel)
            {
                const int x = obj->getHGridX();
                const int y = obj->getHGridY();
                const int z = obj->getHGridZ();
                
                hGridFindContactsWithinTargetCell(x, y, z, level);
                hGridFindContactsWithTargetCell(x + 1, y - 1, z, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y, z, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y + 1, z, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y - 1, z + 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y, z + 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y + 1, z + 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y - 1, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x + 1, y + 1, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x, y + 1, z, level, obj);
                hGridFindContactsWithTargetCell(x, y, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x, y + 1, z - 1, level, obj);
                hGridFindContactsWithTargetCell(x, y + 1, z + 1, level, obj);
            }
            else
            {
                const Mdouble inv_size = hgrid->getInvCellSize(level);
                const int xs = static_cast<int>(std::floor(
                        (obj->getPosition().X - obj->getMaxInteractionRadius()) * inv_size - 0.5));
                const int xe = static_cast<int>(std::floor(
                        (obj->getPosition().X + obj->getMaxInteractionRadius()) * inv_size + 0.5));
                const int ys = static_cast<int>(std::floor(
                        (obj->getPosition().Y - obj->getMaxInteractionRadius()) * inv_size - 0.5));
                const int ye = static_cast<int>(std::floor(
                        (obj->getPosition().Y + obj->getMaxInteractionRadius()) * inv_size + 0.5));
                const int zs = static_cast<int>(std::floor(
                        (obj->getPosition().Z - obj->getMaxInteractionRadius()) * inv_size - 0.5));
                const int ze = static_cast<int>(std::floor(
                        (obj->getPosition().Z + obj->getMaxInteractionRadius()) * inv_size + 0.5));
                for (int x = xs; x <= xe; ++x)
                {
                    for (int y = ys; y <= ye; ++y)
                    {
                        for (int z = zs; z <= ze; ++z)
                        {
                            hGridFindContactsWithTargetCell(x, y, z, level, obj);
                        }
                    }
                }
            }
        }
    }
}

/*!
 * \param[in] obj   A pointer to the BaseParticle that must be updated.
 * \details         Updates the HGrid positions (x, y and z) of the given BaseParticle.
 */
void Mercury3D::hGridUpdateParticle(BaseParticle* obj)
{
    HGrid* const hGrid = getHGrid();
    if (hGrid)
    {
        const unsigned int l = obj->getHGridLevel();
        const Mdouble inv_size = hGrid->getInvCellSize(l);
        
        int x = static_cast<int>(std::floor(obj->getPosition().X * inv_size));
        int y = static_cast<int>(std::floor(obj->getPosition().Y * inv_size));
        int z = static_cast<int>(std::floor(obj->getPosition().Z * inv_size));

#ifdef CONTACT_LIST_HGRID
        if(obj->getHGridX() != x || obj->getHGridY() != y || obj->getHGridZ() != z)
        {
            int bucket = hGrid->computeHashBucketIndex(x, y, z, l);

            //First the object has to be removed
            hGridRemoveParticle(obj);

            //Also remove all contact associated with it
            getPossibleContactList().remove_ParticlePosibleContacts(obj);

            //And now reinserted
            obj->setHGridNextObject(hGrid->getFirstBaseParticleInBucket(bucket));
            obj->setHGridPrevObject(nullptr);
            if(hGrid->getFirstBaseParticleInBucket(bucket))
            {
                hGrid->getFirstBaseParticleInBucket(bucket)->setHGridPrevObject(obj);
            }
            hGrid->setFirstBaseParticleInBucket(bucket,obj);

            obj->setHGridX(x);
            obj->setHGridY(y);
            obj->setHGridZ(z);
            InsertObjAgainstGrid(obj);
        }
#else
        const unsigned int bucket = hGrid->computeHashBucketIndex(x, y, z, l);
        
        // this needs to be defined as #pragma omp critical if MercuryBase::hGridActionsBeforeTimeStep is parallelised; however, parallelising it make the code slower, not faster.
        {
            obj->setHGridNextObject(hGrid->getFirstBaseParticleInBucket(bucket));
            obj->setHGridPrevObject(nullptr);
            if (hGrid->getFirstBaseParticleInBucket(bucket)) {
                hGrid->getFirstBaseParticleInBucket(bucket)->setHGridPrevObject(obj);
            }
            hGrid->setFirstBaseParticleInBucket(bucket, obj);
        }
        
        obj->setHGridX(x);
        obj->setHGridY(y);
        obj->setHGridZ(z);
#endif
    }
}

/*!
 * \param[in] obj   A pointer to the BaseParticle that needs to be removed.
 * \details         Removes the given BaseParticle from the HGrid.
 */
void Mercury3D::hGridRemoveParticle(BaseParticle* obj)
{
    HGrid* const hGrid = getHGrid();
    if (hGrid)
    {
        const unsigned int bucket = hGrid->computeHashBucketIndex(obj->getHGridCell());
        if (obj->getHGridPrevObject())
        {
            obj->getHGridPrevObject()->setHGridNextObject(obj->getHGridNextObject());
        }
        else
        {
            if (hGrid->getFirstBaseParticleInBucket(bucket) == obj)
            {
                hGrid->setFirstBaseParticleInBucket(bucket, obj->getHGridNextObject());
            }
        }
        
        if (obj->getHGridNextObject())
        {
            obj->getHGridNextObject()->setHGridPrevObject(obj->getHGridPrevObject());
        }
    }
}

/*!
 * \param[in] x     The coordinate of the target cell in x direction.
 * \param[in] y     The coordinate of the target cell in y direction.
 * \param[in] z     The coordinate of the target cell in z direction.
 * \param[in] l     The level of the HGrid of the target cell.
 * \param[in] obj   A pointer to the BaseParticle which is checked for contacts.
 * \details Tests if there are any collisions between given BaseParticle and 
 * particles in cell (x, y, z, l).
 */
bool Mercury3D::hGridHasContactsInTargetCell(int x, int y, int z, unsigned int l, const BaseParticle* obj) const
{
    // Loop through all objects in the bucket to find nearby objects
    const unsigned int bucket = getHGrid()->computeHashBucketIndex(x, y, z, l);
    
    const BaseParticle* p = getHGrid()->getFirstBaseParticleInBucket(bucket);
    while (p != nullptr)
    {
        if (p->getHGridCell().equals(x, y, z, l))
        {
            if (areInContact(obj, p))
            {
                return true;
            }
        }
        //std::cout << "HERE!" << std::endl;
        p = p->getHGridNextObject();
    }
    return false;
}

/*!
 * \param[in] obj   A pointer to the BaseParticle that is tested for contacts.
 * \details         Tests if there are any collisions between the given 
 *                  BaseParticle and all other particles in the HGrid. Do this by
 *                  going through all levels, find if there is a collision in any
 *                  of the levels in any cell of the HGrid.
 */
///
bool Mercury3D::hGridHasParticleContacts(const BaseParticle* obj)
{
    if (getHGrid() == nullptr || getHGrid()->getNeedsRebuilding())
    {
        logger(INFO, "HGrid needs rebuilding for \"bool Mercury3D::hGridHasParticleContacts(BaseParticle *obj)\"");
        hGridRebuild();
    }
    
    int occupiedLevelsMask = getHGrid()->getOccupiedLevelsMask();
    
    for (unsigned int level = 0; level < getHGrid()->getNumberOfLevels(); occupiedLevelsMask >>= 1, level++)
    {
        // If no objects in rest of grid, stop now
        if (occupiedLevelsMask == 0)
        {
            logger(VERBOSE, "Level % and higher levels are empty", level);
            break;
        }
        
        // If no objects at this level, go on to the next level
        if ((occupiedLevelsMask & 1) == 0)
        {
            logger(VERBOSE, "Level % is empty", level);
            continue;
        }
        
        const Mdouble inv_size = getHGrid()->getInvCellSize(level);
        const int xs = static_cast<int>(std::floor(
                (obj->getPosition().X - obj->getMaxInteractionRadius()) * inv_size - 0.5));
        const int xe = static_cast<int>(std::floor(
                (obj->getPosition().X + obj->getMaxInteractionRadius()) * inv_size + 0.5));
        const int ys = static_cast<int>(std::floor(
                (obj->getPosition().Y - obj->getMaxInteractionRadius()) * inv_size - 0.5));
        const int ye = static_cast<int>(std::floor(
                (obj->getPosition().Y + obj->getMaxInteractionRadius()) * inv_size + 0.5));
        const int zs = static_cast<int>(std::floor(
                (obj->getPosition().Z - obj->getMaxInteractionRadius()) * inv_size - 0.5));
        const int ze = static_cast<int>(std::floor(
                (obj->getPosition().Z + obj->getMaxInteractionRadius()) * inv_size + 0.5));
        
        logger(VERBOSE, "Level = % grid cells [%,%] x [%,%] x [%,%]", level, xs, xe, ys, ye, zs, ze);
        for (int x = xs; x <= xe; ++x)
        {
            for (int y = ys; y <= ye; ++y)
            {
                for (int z = zs; z <= ze; ++z)
                {
                    if (hGridHasContactsInTargetCell(x, y, z, level, obj))
                    {
                        return true;
                    }
                }
            }
        }
    } //end for level
    
    return false;
}

std::vector<BaseParticle*> Mercury3D::hGridFindParticleContacts(const BaseParticle* obj)
{
    if (getHGrid() == nullptr || getHGrid()->getNeedsRebuilding())
    {
        logger(INFO, "HGrid needs rebuilding for \"bool Mercury3D::hGridHasParticleContacts(BaseParticle *obj)\"");
        hGridRebuild();
    }
    
    int occupiedLevelsMask = getHGrid()->getOccupiedLevelsMask();
    
    std::vector<BaseParticle*> particlesInContact;
    
    for (unsigned int level = 0; level < getHGrid()->getNumberOfLevels(); occupiedLevelsMask >>= 1, level++)
    {
        // If no objects in rest of grid, stop now
        if (occupiedLevelsMask == 0)
        {
            logger(VERBOSE, "Level % and higher levels are empty", level);
            break;
        }
        
        // If no objects at this level, go on to the next level
        if ((occupiedLevelsMask & 1) == 0)
        {
            logger(VERBOSE, "Level % is empty", level);
            continue;
        }
        
        const Mdouble inv_size = getHGrid()->getInvCellSize(level);
        const int xs = static_cast<int>(std::floor(
                (obj->getPosition().X - obj->getMaxInteractionRadius()) * inv_size - 0.5));
        const int xe = static_cast<int>(std::floor(
                (obj->getPosition().X + obj->getMaxInteractionRadius()) * inv_size + 0.5));
        const int ys = static_cast<int>(std::floor(
                (obj->getPosition().Y - obj->getMaxInteractionRadius()) * inv_size - 0.5));
        const int ye = static_cast<int>(std::floor(
                (obj->getPosition().Y + obj->getMaxInteractionRadius()) * inv_size + 0.5));
        const int zs = static_cast<int>(std::floor(
                (obj->getPosition().Z - obj->getMaxInteractionRadius()) * inv_size - 0.5));
        const int ze = static_cast<int>(std::floor(
                (obj->getPosition().Z + obj->getMaxInteractionRadius()) * inv_size + 0.5));
        
        logger(VERBOSE, "Level = % grid cells [%,%] x [%,%] x [%,%]", level, xs, xe, ys, ye, zs, ze);
        for (int x = xs; x <= xe; ++x)
        {
            for (int y = ys; y <= ye; ++y)
            {
                for (int z = zs; z <= ze; ++z)
                {
                    // Loop through all objects in the bucket to find nearby objects
                    const unsigned int bucket = getHGrid()->computeHashBucketIndex(x, y, z, level);
                    BaseParticle* p = getHGrid()->getFirstBaseParticleInBucket(bucket);
                    while (p != nullptr)
                    {
                        if (p->getHGridCell().equals(x, y, z, level))
                        {
                            if (areInContact(obj, p))
                            {
                                particlesInContact.push_back(p);
                            }
                        }
                        p = p->getHGridNextObject();
                    }
                }
            }
        }
    } //end for level
    
    return particlesInContact;
}

void Mercury3D::computeWallForces(BaseWall* const w)
{
    
    // if wall is not local, use the non-hGrid version for finding wall contacts
    Vec3D min, max;
    if (w->isLocal(min, max)==false)
    {
        return DPMBase::computeWallForces(w);
    }
    
    //compute forces for all particles that are neither fixed or ghosts
    if (getHGrid() == nullptr || getHGrid()->getNeedsRebuilding())
    {
        logger(INFO, "HGrid needs rebuilding for \"bool Mercury3D::hGridHasParticleContacts(BaseParticle *obj)\"");
        hGridRebuild();
    }
    
    HGrid* const hGrid = getHGrid();
    
    int occupiedLevelsMask = hGrid->getOccupiedLevelsMask();
    
    for (unsigned int level = 0; level < hGrid->getNumberOfLevels(); occupiedLevelsMask >>= 1, level++)
    {
        // If no objects in rest of grid, stop now
        if (occupiedLevelsMask == 0)
        {
            logger(VERBOSE, "Level % and higher levels are empty", level);
            break;
        }
        
        // If no objects at this level, go on to the next level
        if ((occupiedLevelsMask & 1) == 0)
        {
            logger(VERBOSE, "Level % is empty", level);
            continue;
        }
        
        const Mdouble inv_size = hGrid->getInvCellSize(level);
        const int xs = static_cast<int>(std::floor(min.X * inv_size - 0.5));
        const int xe = static_cast<int>(std::floor(max.X * inv_size + 0.5));
        const int ys = static_cast<int>(std::floor(min.Y * inv_size - 0.5));
        const int ye = static_cast<int>(std::floor(max.Y * inv_size + 0.5));
        const int zs = static_cast<int>(std::floor(min.Z * inv_size - 0.5));
        const int ze = static_cast<int>(std::floor(max.Z * inv_size + 0.5));
        //logger(INFO, "Level % grid cells [%,%] x [%,%] x [%,%]", level, xs, xe, ys, ye, zs, ze);
        
        for (int x = xs; x <= xe; ++x)
        {
            for (int y = ys; y <= ye; ++y)
            {
                for (int z = zs; z <= ze; ++z)
                {
                    // Loop through all objects in the bucket to find nearby objects
                    const unsigned int bucket = hGrid->computeHashBucketIndex(x, y, z, level);
                    BaseParticle* p = hGrid->getFirstBaseParticleInBucket(bucket);
                    while (p != nullptr)
                    {
                        if (!p->isFixed() && p->getPeriodicFromParticle() == nullptr &&
                            p->getHGridCell().equals(x, y, z, level))
                        {
                            //logger(INFO, "t % p % Size % level % cells % % %", getNumberOfTimeSteps(), p->getIndex(), hGrid->getCellSize(level), level, x,y,z);
                            computeForcesDueToWalls(p, w);
                            //w->computeForces(p);
                        }
                        p = p->getHGridNextObject();
                    }
                }
            }
        }
    } //end for level
}

#ifdef CONTACT_LIST_HGRID
/*!
 * \param[in] x     The x coordinate of the cell for which possible contacts are requested.
 * \param[in] y     The y coordinate of the cell for which possible contacts are requested.
 * \param[in] z     The z coordinate of the cell for which possible contacts are requested.
 * \param[in] l     The level in the HGrid of the cell for which possible contacts are requested.
 * \param[in] obj   A pointer to the BaseParticle for which possible contacts are requested.
 * \details         Adds the combination of all objects in the cell with identity (x,y,z,l)
 *                  and given BaseParticle to the list of possible contacts.
 */
void Mercury3D::InsertCell(int x, int y, int z, unsigned int l, BaseParticle *obj)
{   
    // Loop through all objects in the bucket to find nearby objects
    unsigned int bucket = getHGrid()->computeHashBucketIndex(x,y,z,l);
    BaseParticle *p = getHGrid()->getFirstBaseParticleInBucket(bucket);

    while (p!=nullptr)
    {   
        if ((p->getHGridX() == x) && (p->getHGridY() == y) && (p->getHGridZ() == z) && (p->getHGridLevel() == l) && (obj!=p))
        {   
            getPossibleContactList().add_PossibleContact(obj,p);
        }
        p = p->getHGridNextObject();
    }
}

/*!
 * \param[in] obj   A pointer to the BaseParticle for which all possible interactions are requested.
 * \details         Add the object to possible interactions for all levels. 
 *                  Do this by first finding the cell of the given BaseParticle.
 *                  Then go through all levels and find all possible interactions,
 *                  which are then added to the list.
 */
void Mercury3D::InsertObjAgainstGrid(BaseParticle *obj)
{   
    Mdouble inv_size;
    int occupiedLevelsMask_ = getHGrid()->getOccupiedLevelsMask();

    inv_size=getHGrid()->getInvCellSize(obj->getHGridLevel());

    double ownXMin = (obj->getHGridX() - 0.5) * getHGrid()->getCellSize(obj->getHGridLevel());
    double ownXMax = (obj->getHGridX() + 1.5) * getHGrid()->getCellSize(obj->getHGridLevel());
    double ownYMin = (obj->getHGridY() - 0.5) * getHGrid()->getCellSize(obj->getHGridLevel());
    double ownYMax = (obj->getHGridY() + 1.5) * getHGrid()->getCellSize(obj->getHGridLevel());
    double ownZMin = (obj->getHGridZ() - 0.5) * getHGrid()->getCellSize(obj->getHGridLevel());
    double ownZMax = (obj->getHGridZ() + 1.5) * getHGrid()->getCellSize(obj->getHGridLevel());

    for (int level = 0; level < getHGrid()->getNumberOfLevels(); occupiedLevelsMask_ >>= 1, level++)
    {   
        // If no objects in rest of grid, stop now
        if (occupiedLevelsMask_ == 0)
        {   
            break;
        }

        // If no objects at this level, go on to the next level
        if ((occupiedLevelsMask_ & 1) == 0)
        {   
            continue;
        }

        // Treat level as a third dimension coordinate
        inv_size = getHGrid()->getInvCellSize(level);

        int xs, xe, ys, ye, zs, ze;
        xs=static_cast<int>(std::floor(ownXMin * inv_size - 0.5));
        xe=static_cast<int>(std::floor(ownXMax * inv_size + 0.5));
        ys=static_cast<int>(std::floor(ownYMin * inv_size - 0.5));
        ye=static_cast<int>(std::floor(ownYMax * inv_size + 0.5));
        zs=static_cast<int>(std::floor(ownZMin * inv_size - 0.5));
        ze=static_cast<int>(std::floor(ownZMax * inv_size + 0.5));

        for(int x=xs; x<=xe; ++x)
        {   
            for(int y=ys; y<=ye; ++y)
            {   
                for(int z=zs; z<=ze; ++z)
                {   
                    InsertCell(x, y, z, level, obj);
                }
            }
        }
    } //end for level
}
#endif
