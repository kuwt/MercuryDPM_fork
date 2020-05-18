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

#include "HGrid.h"
#include "Logger.h"
#include "Particles/BaseParticle.h"

HGrid::HGrid()
{
    needsRebuilding_ = true;
    numberOfBuckets_ = 10;
    cellOverSizeRatio_ = 1.0;
    occupiedLevelsMask_ = 0;
    logger(DEBUG, "HGrid::HGrid() finished");
}

/*!
 * \param[in] num_buckets       The number of buckets that are used by this HGrid.
 * \param[in] cellOverSizeRatio The maximum ratio between the size of the 
 *                              cell over the size of the particle.
 * \param[in] cellSizes         The sizes of the cells we want to set.
 * \details                     Constructor: initialises parameters and allocates 
 *                              space for internal variables. 
 */
HGrid::HGrid(unsigned int num_buckets, double cellOverSizeRatio, std::vector<double>& cellSizes)
{
    needsRebuilding_ = false;
    numberOfBuckets_ = num_buckets;
    cellOverSizeRatio_ = cellOverSizeRatio;
    occupiedLevelsMask_ = 0;
    invCellSizes_ = std::vector<double>(0);
    
    firstBaseParticleInBucket_.resize(numberOfBuckets_, nullptr);
    bucketIsChecked_.resize(numberOfBuckets_, false);
    
    //std::cout<<"Creating HGrid "<<cellSizes.size()<<" levels:"<<std::endl;
    for (double cellSize : cellSizes)
    {
        //std::cout<<"Level="<<i<<" size="<<cellSizes[i]<<std::endl;
        cellSizes_.push_back(cellSize);
        invCellSizes_.push_back(1.0 / cellSize);
    }
    logger(DEBUG, "HGrid::HGrid(unsigned int, double, vector<double>&) constructor finished.");
    /*  std::cout << "HGrid::HGrid(" << num_buckets << ", " << cellOverSizeRatio << ", [";
        for (auto p: cellSizes) std::cout << p << " "; 
        std::cout << "]) finished" << std::endl;*/
}

HGrid::~HGrid()
{
    logger(DEBUG, "HGrid::~HGrid() destructor finished");
}

/*!
 * \param[in] obj A pointer to the BaseParticle we want to add to the HGrid.
 * \details Inserts the given BaseParticle into the HGrid, i.e. it sets up the 
 *          particle grid properties and updates the level information on the grid.
 *          First find which level is big enough to fit the BaseParticle in, then 
 *          add the BaseParticle to that level and set that level as occupied in 
 *          the occupiedLevelsMask_.
 * \bug What happens if the particle is too big for the biggest cell? It just says
 *      that it needs to rebuild the HGrid, but the particle is not inserted and
 *      there seems to be no indication to the rest of the code that it has not 
 *      been inserted. For now giving a warning, since code of users may rely on
 *      it that nothing happens.
 */
void HGrid::insertParticleToHgrid(BaseParticle* obj)
{
    if (!needsRebuilding_)
    {
        // Find lowest level where object fully fits inside cell, taking cellOverSizeRatio_ into account
        Mdouble diameter = obj->getMaxInteractionRadius() * 2.0;
        unsigned int level = 0;
        while (level < (cellSizes_.size() - 1) && cellSizes_[level] <= diameter * cellOverSizeRatio_)
        {
            level++;
        }
        
        //Check if the size of the particle is larger than the required grid
        if (level >= cellSizes_.size())
        {
            logger(WARN, "WARNING: object (id = %, index = %) is larger (d = %, cellOverSizeRatio = %) than largest "
                         "grid cell (%) allows.",
                   obj->getId(), obj->getIndex(), diameter, cellOverSizeRatio_, cellSizes_.back());
            needsRebuilding_ = true;
        }
        
        obj->setHGridLevel(level);
        // indicate level is in use - not levels with no particles no collision detection is performed
        this->occupiedLevelsMask_ |= (1 << level);
    }
    else
    {
        logger(WARN, "WARNING: the HGrid needs to be rebuild before insertParticleToHgrid may be called!");
    }
}

/*!
 * \details Computes a hash from parameters, the result is in range [0, numberOfBuckets_-1].
 * \param[in] x The coordinate of the cell in x direction for which the hash must be computed.
 * \param[in] y The coordinate of the cell in y direction for which the hash must be computed.
 * \param[in] l The level in the HGrid of the cell for which the hash must be computed.
 * \return The hash value for the given cell (x,y,l), which is in the range [0,numberOfBuckets_-1].
 */
unsigned int HGrid::computeHashBucketIndex(int x, int y, unsigned int l) const
{
    const unsigned int h1 = 0x8da6b343u; // Large multiplicative constants;
    const unsigned int h2 = 0xd8163841u; // here arbitrarily chosen primes
    const unsigned int h4 = 0x165667b1u;
    
    unsigned long int n = h1 * x + h2 * y + h4 * l;
    n = n % numberOfBuckets_;
    
    return static_cast<unsigned int>(n);
}


void HGrid::clearBucketIsChecked()
{
    std::fill(bucketIsChecked_.begin(), bucketIsChecked_.end(), false);
}

void HGrid::clearFirstBaseParticleInBucket()
{
    std::fill(firstBaseParticleInBucket_.begin(), firstBaseParticleInBucket_.end(), nullptr);
}

///\todo use logger everywhere
void HGrid::info() const
{
    logger(INFO, "  numberOfBuckets %", numberOfBuckets_);
    logger(INFO, "  cellOverSizeRatio %", cellOverSizeRatio_);
    std::cout << "  cellSizes";
    for (auto p: cellSizes_) std::cout << " " << p;
    std::cout << '\n';
}
