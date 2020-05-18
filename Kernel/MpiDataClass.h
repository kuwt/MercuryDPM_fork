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

#ifndef MPIDATACLASS_H_
#define MPIDATACLASS_H_

#include <Particles/SphericalParticle.h>
#include "Particles/BaseParticle.h"
#include "ParticleHandler.h"

class SuperQuadricParticle;

/*!
 * \class MPIParticle
 * \brief Data class to send a particle over MPI
 */
class MPISphericalParticle
{
public:
    unsigned int id;
    unsigned int indSpecies;
    unsigned int HGridLevel;
    Mdouble radius;
    Vec3D position;
    Vec3D angularVelocity;
    Vec3D velocity;
    Quaternion orientation;
    unsigned communicationComplexity;
    bool isMaser; //TODO 
    bool isFixed;
    void copyDataFromMPIParticleToParticle(BaseParticle* p);
    void copyDataFromParticleToMPIParticle(BaseParticle* p);
    static BaseParticle* newParticle ();
};

class MPISuperQuadric : public MPISphericalParticle
{
public:
    Vec3D axes;
    Mdouble epsilon1;
    Mdouble epsilon2;
    void copyDataFromMPIParticleToParticle(BaseParticle* p);
    void copyDataFromParticleToMPIParticle(BaseParticle* p);
    static BaseParticle* newParticle ();
};

class MPILiquidFilmParticle : public MPISphericalParticle
{
public:
    Mdouble liquidVolume;
    void copyDataFromMPIParticleToParticle(BaseParticle* p);
    void copyDataFromParticleToMPIParticle(BaseParticle* p);
    static BaseParticle* newParticle ();
};

/*!
 * Define what type of particle is used in MPI
 */
class MPIParticle : public MPISphericalParticle {};
//to run simulations with LiquidFilmParticles in parallel, uncomment the line below (and comment the line above).
//class MPIParticle : public MPILiquidFilmParticle {};

/*!
 * \class MPIParticlePosition
 * \brief Data class to send a particle position over MPI
 */
class MPIParticlePosition
{
public:
    unsigned int id;
    Vec3D position;
    Quaternion orientation;
    Mdouble liquidVolume;
};

/*!
 * \class MPIParticleVelocity
 * \brief Data class to send a particle velocity over MPI
 */
class MPIParticleVelocity
{
public:
    Vec3D velocity;
    Vec3D angularVelocity;
};

/*!
 * \class MPIParticleForce
 * \brief Data class to send a particle force over MPI
 */
class MPIParticleForce
{
public:
    Vec3D force;
    Vec3D torque;
};

/*!
 * \class Empty
 * \brief Data class to send an empty class over MPI
 * \details The interaction data class has to be flexible since it should work
 * for all valid interaction types (adhesive etc). These interactions are equiped with
 * different history values such as torsionSpring, glued or liquidbridgeVolume.
 * For reducing the amount of data communicated over MPI,
 * these values that are not required for the interaction are
 * replaced by this empty class
 */
class Empty
{
};

/*!
 * \class MpiID
 * \brief Data class that specifies the location of a particle in a parallel code
 * \details With a processor and an id, a particle has a unique position. This is especially
 * used in parallel periodic boundaries where the specific particle particle needs to know where his
 * ghost is hanging around.
 * previousPosition is the position of the particle in the previous time step. This is used to update the particle
 * the reason the previousPosition of the BaseParticle is _not_ used is that that variable is occasionally used in CGing.
 * mixing the two together is not a good idea
 * TODO this documentation needs a severe update
 * TODO it is neater if I additionally create an MpiPeriodicGhostID
 * TODO when mpi delets a particle that is a ghost particle, the periodic list needs to be flushed as well
 */
class MpiPeriodicParticleIDBase
{
public:
    BaseParticle* particle; //Used to store the pointer of the current particle (ghost or periodic) located on the same domain
    BaseParticle* otherParticle; //Used to store the pointer of the other particle (ghost or periodic) located on the same domain
    //int currentProcessor;
    //int otherProcessor;
    int targetProcessor; //Stores the target processor where the ghost or periodic is located
    //int targetID;
    //int periodicBoundaryIndex;
    //std::vector<int> previousPeriodicComplexity;
    //std::vector<int> currentPeriodicComplexity;
    std::vector<int> periodicComplexity; //Not sure if this is required
    std::vector<int> targetPeriodicComplexity; //ppid uses this to store the ghost periodic complexity when adding
    std::vector<int> realPeriodicComplexity; //gpid uses this to store the realPeriodicComplexity of the real particle
    //std::vector<int> previousRealPeriodicComplexity;
    //Vec3D targetPosition;
};

typedef MpiPeriodicParticleIDBase MpiPeriodicParticleID;
typedef MpiPeriodicParticleIDBase MpiPeriodicGhostParticleID;

/*!
 * \brief Copies data from a BaseParticle to an MPIParticle class and returns this
 */
MPIParticle copyDataFromParticleToMPIParticle(BaseParticle* p);

/*!
 * \brief Copies data from an MPIParticle class to a BaseParticle
 */
void copyDataFromMPIParticleToParticle(MPIParticle* bP, BaseParticle* p, ParticleHandler* particleHandler);

/*!
 * \brief Copies the position from a particle to an MPIParticlePosition class
 */
MPIParticlePosition copyPositionFrom(BaseParticle* particle);

/*!
 * \brief Copies the velocity from a particle to an MPIParticleVelocity class
 */
MPIParticleVelocity copyVelocityFrom(BaseParticle* particles);

/**
 * Sums the values over all processors using MPI_reduce
 */
Vec3D getMPISum(Vec3D& val);

/**
 * Sums the values over all processors using MPI_reduce
 */
double getMPISum(double val);

#endif  /* MPIDATACLASS_H_ */
