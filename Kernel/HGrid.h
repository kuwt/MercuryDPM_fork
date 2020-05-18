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

#ifndef HGRID_H
#define HGRID_H

#include <vector>
#include "Math/ExtendedMath.h"
#include "Particles/HGridCell.h"

class BaseParticle;

/*!
 * \brief In the HGrid class, here all information about the HGrid is stored.
 * \details In particular, the hashing the grid is done in this class, the cell 
 *          sizes of the different levels are stored, it is stored whether or not
 *          there are BaseParticle on each level and there is a flag to see if the
 *          HGrid needs to be rebuilt.
 */
class HGrid
{
public:
    /*!
     * \brief Default constructor, it sets the parameters to some sensible defaults.
     */
    HGrid();
    
    /*!
     * \brief Constructor: initialises parameters and allocates space for internal variables
     */
    HGrid(unsigned int num_buckets, double cellOverSizeRatio, std::vector<double>& cellSizes);
    
    /*!
     * \brief Destructor.
     */
    ~HGrid();
    
    /*!
     * \brief Inserts the given BaseParticle in to the HGrid.
     */
    void insertParticleToHgrid(BaseParticle* obj);
    
    /*!
     * \brief Computes hash bucket index in range [0, NUM_BUCKETS-1] for a 3D domain.
     * \details Computes a hash from parameters, the result is in range [0, numberOfBuckets_-1].
     * Inline for performance reasons: this method is called VERY often.
     * \param[in] x The coordinate of the cell in x direction for which the hash must be computed.
     * \param[in] y The coordinate of the cell in y direction for which the hash must be computed.
     * \param[in] z The coordinate of the cell in z direction for which the hash must be computed.
     * \param[in] l The level in the HGrid of the cell for which the hash must be computed.
     * \return The hash value for the given cell (x,y,z,l), which is in the range [0,numberOfBuckets_-1].
     * \todo consider moving to HGridCell, which might give better performance.
     */
    inline unsigned int computeHashBucketIndex(int x, int y, int z, unsigned int l) const
    {
        static const unsigned int h1 = 0x8da6b343u; // Large multiplicative constants;
        static const unsigned int h2 = 0xd8163841u; // here arbitrarily chosen primes
        static const unsigned int h3 = 0xcb1ab31fu;
        static const unsigned int h4 = 0x165667b1u;
        
        return (h1 * x + h2 * y + h3 * z + h4 * l) % numberOfBuckets_;
    }
    
    /*!
     * \brief Computes hash bucket index in range [0, NUM_BUCKETS-1] for a 3D domain.
     */
    inline unsigned int computeHashBucketIndex(HGridCell hGridCell) const
    {
        return computeHashBucketIndex(hGridCell.getHGridX(),
                                      hGridCell.getHGridY(),
                                      hGridCell.getHGridZ(),
                                      hGridCell.getHGridLevel());
    }
    
    /*!
     * \brief Computes hash bucket index in range [0, NUM_BUCKETS-1] for a 2D domain.
     */
    unsigned int computeHashBucketIndex(int x, int y, unsigned int l) const;
    
    /*!
     * \brief Sets all buckets to not-checked.
     */
    void clearBucketIsChecked();
    
    /*!
     * \brief For all buckets, it removes the pointer to the first BaseParticle in it, practically emptying the buckets.
     */
    void clearFirstBaseParticleInBucket();
    
    /*!
     * \brief Sets the first particle in bucket i to be the given BaseParticle.
     * \param[in] i The ordinal number of the bucket we want to set the first BaseParticle for.
     * \param[in] p A pointer to the BaseParticle we want to place in the given bucket.
     */
    void setFirstBaseParticleInBucket(unsigned int i, BaseParticle* p)
    { firstBaseParticleInBucket_[i] = p; }
    
    /*!
     * \brief Sets that the bucket with the given index is checked to true.
     * \param[in] i The ordinal number of the bucket we want to mark as checked.
     */
    void setBucketIsChecked(unsigned int i)
    { bucketIsChecked_[i] = true; }
    
    /*!
     * \brief Gets whether or not the bucket with index i is checked.
     * \param[in] i The ordinal number of the bucket we want to know for whether or not it has been checked.
     * \return A boolean which is true if the bucket is checked and false otherwise.
     */
    bool getBucketIsChecked(unsigned int i) const
    { return bucketIsChecked_[i]; }
    
    /*!
     * \brief Gets the maximum ratio of the cell to a particle it contains.
     * \return The ratio between the size of the smallest cell and the smallest BaseParticle.
     */
    Mdouble getCellOverSizeRatio() const
    { return cellOverSizeRatio_; }
    
    /*!
     * \brief Gets the size of the cells at the given level.
     * \param[in] i The level we want to know the cell size of.
     * \return The size of the cells at the given level.
     */
    double getCellSize(unsigned int i) const
    { return cellSizes_[i]; }
    
    /*!
     * \brief Gets the sizes of the cells at all levels as a vector.
     * \return A vector with the sizes of the cells of different levels.
     */
    const std::vector<double>& getCellSizes() const
    { return cellSizes_; }
    
    /*!
     * \brief Gets the first BaseParticle in the given bucket, const version.
     * \param[in] i The ordinal number of the bucket for which we want to get the first particle of.
     * \return A pointer to the (constant) BaseParticle which is the first Baseparticle in the given bucket.
     */
    const BaseParticle* getFirstBaseParticleInBucket(unsigned int i) const
    { return firstBaseParticleInBucket_[i]; }
    
    /*!
     * \brief Gets the first BaseParticle in the given bucket.
     * \param[in] i The ordinal number of the bucket for which we want to get the first particle of.
     * \return A pointer to the (constant) BaseParticle which is the first Baseparticle in the given bucket.
     */
    BaseParticle* getFirstBaseParticleInBucket(unsigned int i)
    { return firstBaseParticleInBucket_[i]; }
    
    /*!
     * \brief Gets 1/cellSize for the cells on level i.
     * \param[in] i The level we want to know the inverse cell size of.
     * \return The inverse size, i.e. 1/size, of the cells at the given level.
     */
    double getInvCellSize(unsigned int i) const
    { return invCellSizes_[i]; }
    
    /*!
     * \brief Gets all the inverse cell sizes (1/cellSize) for all levels as a vector.
     * \return A vector with the inverse sizes (1/size) of the cells of different levels.
     */
    const std::vector<double>& getInvCellSizes() const
    { return invCellSizes_; }
    
    /*!
     * \brief Gets whether or not the grid needs to be rebuilt before something else is done with it.
     * \return A boolean which indicates whether or not the HGrid needs rebuilding.
     */
    bool getNeedsRebuilding() const
    { return needsRebuilding_; }
    
    /*!
     * \brief Gets the number of buckets of this HGrid.
     * \return The number of buckets in this HGrid.
     */
    unsigned int getNumberOfBuckets() const
    { return numberOfBuckets_; }
    
    /*!
     * \brief Gets the number of levels of this HGrid.
     * \return The number of levels in this HGrid.
     */
    unsigned long getNumberOfLevels() const
    { return cellSizes_.size(); }
    
    /*!
     * \brief Gets the integer that represents which levels are occupied.
     * \return The integer that represents the bit-vector that indicates which levels have at least one particle.
     */
    int getOccupiedLevelsMask() const
    { return occupiedLevelsMask_; }
    
    /*!
     * \brief Displays the member variables of the hGrid object.
     * This function is intended for debugging the hGrid, therefore the 
     * variables are displayed next to the variable names instead of putting it 
     * in a user-friendly format. 
     */
    void info() const;

private:
    /*!
     * \brief Flag sets if the HGrid needs to be rebuilt.
     */
    bool needsRebuilding_;
    
    /*!
     * \brief The number of buckets in the current HGrid.
     * \details The number of buckets is the number of possible "results" of the
     *          hash function for the grid.
     */
    unsigned int numberOfBuckets_;
    
    /*!
     * \brief The maximum ratio between the size of the cell and the size of a particle it contains.
     */
    Mdouble cellOverSizeRatio_;
    
    /*!
     * \brief Marks if there are particles at certain levels.
     * \details The l-th bit of occupiedLevelsMask_ is 1 if level l is contains 
     *          particles (Implies max 32 hgrid levels) and 0 if it contains none.
     */
    int occupiedLevelsMask_;
    
    /*!
     * \brief The sizes of the cells in the different grids.
     * \details The sizes of the cells in the different grids are saved in a 
     *          vector of double. The smaller the index in the vector, the smaller
     *          the cells.
     */
    std::vector<double> cellSizes_;
    
    /*!
     * \brief The inverse sizes of the cells in the different grids, where the inverse is defined as 1/cellSizes_.
     */
    std::vector<double> invCellSizes_;
    
    /*!
     * \brief Stores a pointer to first element in hash bucket b.
     * \details The pointer to the first element in a certain bucket, initially
     *          a nullptr, is the pointer to the first BaseParticle in the first
     *          cell in the bucket.
     */
    std::vector<BaseParticle*> firstBaseParticleInBucket_;
    
    /*!
     * \brief BucketIsChecked stores if hash bucket b is checked already; initially all false.
     */
    std::vector<bool> bucketIsChecked_;
};

#endif
