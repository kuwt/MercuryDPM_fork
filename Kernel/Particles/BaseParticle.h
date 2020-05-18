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

#ifndef BASEPARTICLE_H
#define BASEPARTICLE_H

#include <vector>
#include "Math/Vector.h"
#include "BaseInteractable.h"
#include "HGridCell.h"
#include "Species/ParticleSpecies.h"

class ParticleHandler;

class SpeciesHandler;

class HGrid;

class BaseInteraction;

#ifdef CONTACT_LIST_HGRID
class PossibleContact;
#endif

/*!
 * \class BaseParticle
 * \brief
 * \details Since r3648, BaseParticle is an abstract class. Use
 * SphericalParticle for a 'basic' particle.
 */
class BaseParticle : public BaseInteractable
{
public:
    /*!
     * \brief Basic Particle constructor, creates an Particle at (0,0,0) with radius, mass and inertia equal to 1
     */
    BaseParticle();
    
    /*!
     * \brief Particle copy constructor, which accepts as input a reference to a Particle. It creates a copy of this Particle and all it's information. Usually it is better to use the copy() function for polymorfism.
     */
    BaseParticle(const BaseParticle& p);
    
    explicit BaseParticle(const ParticleSpecies* s);
    
    /*!
     * \brief Particle destructor, needs to be implemented and checked if it removes tangential spring information
     */
    ~BaseParticle() override;
    
    /*!
     * \brief Particle copy method. It calls to copy constructor of this Particle, useful for polymorfism
     */
    virtual BaseParticle* copy() const = 0;

    /*!
     * \brief Get Particle volume function, which required a reference to the Species vector. It returns the volume of the Particle.
     */
    virtual Mdouble getVolume() const;
    
    /*!
     * \brief Fix Particle function. It fixes a Particle by setting its inverse mass and inertia and velocities to zero
     */
    void fixParticle();
    
    /*!
     * \brief Is fixed Particle function. It returns whether a Particle is fixed or not, by checking its inverse Mass.
     * \details Checks whether a BaseParticle is fixed or not, by checking its inverse Mass.
     * \return TRUE if particle is fixed, i.e. if the inverse mass (invMass_) is 0.
     */
    bool isFixed() const override
    { return (invMass_ == 0.0); }
    
    /*!
     * \brief Indicates if this particle is a ghost in the MPI domain
     */
    bool isMPIParticle() const;
    
    /*!
     * \brief Flags the mpi particle status
     */
    void setMPIParticle(bool flag);
    
    /*!
     * \brief Indicates if the particle is in the communication zone of the mpi domain
     */
    bool isInMPIDomain();
    
    /*!
     * \brief Flags the status of the particle if wether it is in the communication zone or not
     */
    void setInMPIDomain(bool flag);
    
    /*!
     * \brief Indicates if the particle is in the periodic boundary communication zone
     */
    bool isInPeriodicDomain() const;
    
    /*!
     * \brief Flags the status of the particle wether it is in the periodic communication zone or not
     */
    void setInPeriodicDomain(bool flag);
    
    /*!
     * \brief Indicates if this particle is a ghost in the periodic boundary
     */
    bool isPeriodicGhostParticle() const;
    
    /*!
     * \brief Flags the status of the particle to be a ghost in periodic boundary or not
     */
    void setPeriodicGhostParticle(bool flag);
    
    /*!
     * \brief Indicates if this particle belongs to the maser boundary
     */
    bool isMaserParticle() const;
    
    /*!
     * \brief Flags the status of the particle if it belongs to the maser boundary or not
     */
    void setMaserParticle(bool flag);
    
    /*!
     * \brief Set the communication complexity of the particle
     */
    void setCommunicationComplexity(unsigned complexity);
    
    /*!
     * \brief Obtains the communication complexity of the particle
     */
    unsigned getCommunicationComplexity();
    
    /*!
     * \brief Set the periodic communication complexity of the particle
     */
    void setPeriodicComplexity(std::vector<int> complexity);
    
    /*!
     * \brief Set the periodic communication complexity of the particle
     */
    void setPeriodicComplexity(int index, int value);
    
    /*!
     * \brief Obtains the periodic communication complexity of the particle
     */
    const std::vector<int>& getPeriodicComplexity();
    
    /*!
     * \brief Set the previous periodic communication complexity of the paritcle
     */
    void setPreviousPeriodicComplexity(std::vector<int> complexity);
    
    /*!
     * \brief Sets the previous periodic communication complexity of the particle
     */
    const std::vector<int>& getPreviousPeriodicComplexity() const;
    
    /*!
     * \brief Gets the periodic communication complexity of a certain boundary
     */
    int getPeriodicComplexity(int index);
    
    //void setCurrentPeriodicPosition(Vec3D position);
    
    //Vec3D getCurrentPeriodicPosition();   
    
    /*!
     * \brief Unfix Particle function, which required a reference to the Species vector. It unfixes a Particle by computing the Particles mass and inertia
     */
    void unfix();
    
    /*!
     * \brief Particle read function, which accepts an std::istream as input.
     */
    void read(std::istream& is) override;
    
    /*!
     * Should NOT BE USED by any user, only used to read old restart files!
     * Is expected to be obsolete by Mercury 2.0. Please use BaseParticle::read() instead.
     */
    virtual void oldRead(std::istream& is);
    
    /*!
     * \brief Particle print function, which accepts an std::ostream as input.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;
    
    /*!
     * \brief Sets some user-defined information about this object (by default, species ID). 
     */
    virtual void setInfo(Mdouble info);
    
    /*!
     * \brief Returns some user-defined information about this object (by default, species ID).
     */
    virtual Mdouble getInfo() const;
    
    /*!
     * \brief Adds particle's HGrid level and cell coordinates to an ostream
     */
    void printHGrid(std::ostream& os) const;
    
    /*!
     * \brief Returns particle's HGrid level
     */
    unsigned int getHGridLevel() const
    { return hGridCell.getHGridLevel(); }
    
    /*!
     * \brief Returns pointer to next object in particle's HGrid level & cell
     * \details Returns the next object in the particle's HGrid cell
     * \return pointer to the next object in the particle's HGrid cell
     */
    BaseParticle* getHGridNextObject() const
    { return hGridNextObject_; }
    
    /*!
     * \brief Returns pointer to previous object in particle's HGrid level & cell
     * \details Returns the previous object in the particle's HGrid cell
     * \return pointer to the previous object in the particle's HGrid cell
     */
    BaseParticle* getHGridPrevObject() const
    { return hGridPrevObject_; }

#ifdef CONTACT_LIST_HGRID
    
    /*!
     * \brief 
     */
    PossibleContact* getFirstPossibleContact() const;
#endif
    
    /*!
     * \brief Returns particle's HGrid cell X-coordinate
     * \return the particle's HGrid cell's X-coordinate
     */
    int getHGridX() const
    { return hGridCell.getHGridX(); }
    
    /*!
     * \brief Returns particle's HGrid cell Y-coordinate
     * \return the particle's HGrid cell's Y-coordinate
     */
    int getHGridY() const
    { return hGridCell.getHGridY(); }
    
    /*!
     * \brief Returns particle's HGrid cell Z-coordinate
     * \return the particle's HGrid cell's Z-coordinate
     */
    int getHGridZ() const
    { return hGridCell.getHGridZ(); }
    
    /*!
     * \brief Returns the inverse of the particle's inertia tensor
     * \return the inverse of the particle's inertia tensor
     */
    MatrixSymmetric3D getInvInertia() const
    { return invInertia_; }
    
    /*!
     * \brief Returns the inverse of the particle's mass
     * \return the inverse of the particle's mass
     */
    Mdouble getInvMass() const override
    { return invMass_; }

    Mdouble getCurvature(const Vec3D& labFixedCoordinates) const override
    { return 1.0/radius_; }

    /*!
     * \brief Calculates the particle's translational kinetic energy
     * \return the particle's translational kinetic energy
     */
    Mdouble getKineticEnergy() const;
    
    /*!
     * \brief Calculates the particle's rotational kinetic energy
     * \return the particle's rotational kinetic energy
     */
    Mdouble getRotationalEnergy() const;
    
    /*!
     * \brief Calculates the particle's gravitational energy
     * \details Gravitational energy is the potential energy stored in teh particles position due to the gravity field. This is a relative term, so we need to define what zero energy means: The gravitational energy of a particle is zero when its center of mass is at the origin.
     * \return the particle's gravitational energy
     */
    Mdouble getGravitationalEnergy() const;
    
    /*!
     * \brief Returns the particle's mass
     * \return the mass of the particle
     */
    Mdouble getMass() const
    { return 1.0 / invMass_; }
    
    Vec3D getMomentum() const
    { return getVelocity() / invMass_; }

    MatrixSymmetric3D getInertia() const
    { return invInertia_.inverse(); }
    
    Vec3D getAngularMomentum() const;

    /*!
     * \brief Returns the 'original' particle this one's a periodic copy of
     * \details Returns a pointer to the 'original' particle if the current one is a 'periodic copy' used for a periodic boundary implementation.
     \return pointer to original particle
     */
    BaseParticle* getPeriodicFromParticle() const
    { return periodicFromParticle_; }
    
    /*!
     * \brief Returns the particle's radius
     * \return the particle's radius
     */
    Mdouble getRadius() const
    { return radius_; }
    
    /*!
     * \brief Returns the particle's interaction radius, which might be different
     * from radius_ (e.g., when dealing with wet particles)
     * \details Calculates the interaction radius of the particle (when it comes to
     * interaction with other particles), including the effect of a possible additional
     * 'interaction distance' besides the 'normal' radius. The interaction radius
     * differs from the radius_, for example, when dealing with wet particles (i.e.
     * particles with an additional liquid layer, which is dealt with in the particle's
     * species).
     * \return the particle's interaction radius for particle-particle interaction
     */
    Mdouble getMaxInteractionRadius() const {
        return getRadius() + getSpecies()->getMaxInteractionDistance() * 0.5;
    }
    
    /**
     * \brief Returns the interactionDistance_ of the mixed species of this particle and the particle or wall i
     */
    Mdouble getInteractionDistance(const BaseInteractable* i) const
    {
        //const auto mixedSpecies = getSpecies()->getHandler()->getMixedObject(getSpecies(),particle->getSpecies());
        //return mixedSpecies->getInteractionDistance();
        return getSpecies()->getMixedSpecies(i->getSpecies())->getInteractionDistance();
    }
    
    /**
     * \brief returns the sum of the radii plus the interactionDistance
     */
    Mdouble getSumOfInteractionRadii(const BaseParticle* particle) const {
        return getRadius() + particle->getRadius() + getInteractionDistance((const BaseInteractable*)particle);
    }
    
    /**
     * \brief returns the radius plus the interactionDistance
     */
    Mdouble getWallInteractionRadius(const BaseWall* wall) const {
        return getRadius() + getInteractionDistance((const BaseInteractable*)wall);
    }
    
    /*!
     * \brief Returns the particle's displacement relative to the previous time step
     * \details Returns the particle's displacement_, which is the difference between the current particle's position and its position in the previous time step.
     * \return (reference to) the particle displacement vector
     */
    const Vec3D& getDisplacement() const
    { return displacement_; }
    
    /*!
     * \brief Returns the particle's position in the previous time step
     * \details Returns the particle's position in the previous time step.
     * \return (reference to) the previous position of the particle
     */
    const Vec3D& getPreviousPosition() const
    { return previousPosition_; }
    
    /*!
     * \brief 
     * \todo see .cc file. 
     * \TWH
     */
    const Vec3D getDisplacement2(Mdouble xmin, Mdouble xmax, Mdouble ymin, Mdouble ymax, Mdouble zmin, Mdouble zmax,
                                 Mdouble t) const;
    
    /*!
     * \brief 
     */
    virtual void setInertia();
    
    /*!
     * \brief Sets the particle's inertia_ (and adjusts invInertia_ accordingly)
     */
    void setInertia(MatrixSymmetric3D inertia);
    
    /*!
     * \brief Sets the particle's inertia_ (and adjusts invInertia_ accordingly)
     */
    void setInverseInertia(MatrixSymmetric3D inverseInertia);
    
    /*!
     * \brief Sets the particle's inertia_ to 'infinite' (1e20) and its invInertia_
     * to 0.
     */
    void setInfiniteInertia();
    
    /*!
     * \brief Assigns the pointer to the 'original' particle this one's a 
     * periodic copy of (used in periodic boundary condition implementations).
     * \param[in] p  pointer to the 'original' particle this one is a periodic copy of.
     *
     */
    void setPeriodicFromParticle(BaseParticle* p)
    { periodicFromParticle_ = p; };
    
    /*!
     * \brief Sets the particle's HGrid cell X-coordinate
     * \details Set the x-index of the particle's hGrid cell position
     * \param[in] x  x-index of particle's HGrid cell
     */
    void setHGridX(const int x)
    { hGridCell.setHGridX(x); }
    
    /*!
     * \brief Sets the particle's HGrid cell Y-coordinate
     * \details Set the y-index of the particle's hGrid cell position
     * \param[in] y  y-index of particle's HGrid cell
     */
    void setHGridY(const int y)
    { hGridCell.setHGridY(y); }
    
    /*!
     * \brief Sets the particle's HGrid cell Z-coordinate
     * \details Set the y-index of the particle's hGrid cell position
     * \param[in] z  z-index of particle's HGrid cell
     */
    void setHGridZ(const int z)
    { hGridCell.setHGridZ(z); }
    
    /*!
     * \brief Sets the particle's HGrid level
     * \param[in] level  the particle's HGrid level
     */
    void setHGridLevel(const unsigned int level)
    { hGridCell.setHGridLevel(level); }
    
    /*!
     * \brief Sets the pointer to the next object in the particle's HGrid cell
     * & level.
     * \param[in] p     pointer to the next object
     */
    void setHGridNextObject(BaseParticle* p)
    { hGridNextObject_ = p; }
    
    /*!
     * \brief Sets the pointer to the previous object in the particle's HGrid cell
     * & level.
     * \param[in] p     pointer to the previous object
     */
    void setHGridPrevObject(BaseParticle* p)
    { hGridPrevObject_ = p; }

#ifdef CONTACT_LIST_HGRID
    /*!
     * \brief 
     */
    void setFirstPossibleContact(PossibleContact* PC);
#endif
    
    /*!
     * \brief Sets the particle's radius_ (and adjusts the mass_ accordingly, 
     * based on the particle's species)
     */
    virtual void setRadius(Mdouble radius);
    
    /*
     * The below 6 functions had to be declared virtual here in order to allow access of these functions in ParticleHandler class in the function writeVTK.
     */
    
    /*!
     * \brief Only ustilised in case of superquadric particles. Had to create a virtual function to allow function access in writeVTK function in the particle handler.
     */
    virtual Vec3D getAxes() const;
    
    /*!
     * \brief Only ustilised in case of superquadric particles. Had to create a virtual function to allow function access in writeVTK function in the particle handler.
     */
    virtual Mdouble getExponentEps1() const;
    
    /*!
     * \brief Only ustilised in case of superquadric particles. Had to create a virtual function to allow function access in writeVTK function in the particle handler.
     */
    virtual Mdouble getExponentEps2() const;
    
    
    /*!
     * \brief Only ustilised in case of superquadric particles.
     */
    virtual void setAxes(const Vec3D& axes){ }
    
    /*!
     * \brief Only ustilised in case of superquadric particles.
     */
    virtual void setExponents(const Mdouble& eps1, const Mdouble& eps2){}
    
    //
    /*!
     * \brief Sets the particle's mass
     * \deprecated Please do not set the mass yourself, but use 
     * ParticleSpecies->computeMass instead. That makes sure 
     */
    MERCURY_DEPRECATED
    void setMass(Mdouble mass);
    
    /*!
     * \brief Sets the particle's mass
     * This function should not be used, but is necessary to extend the CG 
     * toolbox to non-spherical particles.
     */
    void setMassForP3Statistics(Mdouble mass);
    
    /*!
     * \brief Sets the particle's displacement (= difference between current 
     * position and that of the previous time step)
     */
    void setDisplacement(const Vec3D& disp);
    
    /*!
     * \brief Sets the particle's position in the previous time step
     */
    void setPreviousPosition(const Vec3D& pos);
    
    /*!
     * \brief Adds a vector to the particle's previousPosition_
     */
    void movePrevious(const Vec3D& posMove);
    
    /*!
     * \brief Increases the particle's velocity_ by the given vector
     */
    void accelerate(const Vec3D& vel);
    
    /*!
     * \brief Increases the particle's angularVelocity_ by the given vector
     */
    void angularAccelerate(const Vec3D& angVel);
    
    /*!
     * \brief Adds a vector to the particle's displacement_
     */
    void addDisplacement(const Vec3D& addDisp);
    
    /*!
     * \brief Sets the pointer to the particle's ParticleHandler
     */
    void setHandler(ParticleHandler* handler);
    
    /*!
     * \brief Returns pointer to the particle's ParticleHandler
     */
    ParticleHandler* getHandler() const;
    
    /*!
     * \brief Checks if particle is in interaction with given particle P, and if 
     * so, returns vector of pointer to the associated BaseInteraction object (else returns empty vector).
     */
    BaseInteraction* getInteractionWith(BaseParticle* P, unsigned timeStamp, InteractionHandler* interactionHandler) override;
    
    /*!
     * \brief Get whether or not this particle is in contact with the given particle.
     */
    virtual bool isInContactWith(const BaseParticle* P) const;
    
    /*!
     * \brief First step of <a href="http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet">Velocity Verlet</a> 
     * integration.
     * 
     */
    void integrateBeforeForceComputation(double time, double timeStep);
    
    /*!
     * \brief Second step of <a href="http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet">Velocity Verlet</a> 
     * integration.
     */
    void integrateAfterForceComputation(double time, double timeStep);
    
    /*!
     * \brief Returns the particle's dimensions (either 2 or 3).
     */
    unsigned int getParticleDimensions() const;
    
    /*!
     * \deprecated Please use setSpecies(const ParticleSpecies*) instead.
     */
    MERCURY_DEPRECATED
    /// \todo MX: this index is used in the MPI transmission. This should be "undeprecated"
    void setIndSpecies(unsigned int indSpecies) override;
    
    /*!
     * In addition to the functionality of BaseInteractable::setSpecies, this function sets the pointer to the
     * particleHandler, which is needed to retrieve species information.
     */
    ///\todo TW: this function should also check if the particle is the correct particle for the species type
    void setSpecies(const ParticleSpecies* species);
    
    virtual unsigned getNumberOfFieldsVTK() const;
    
    virtual std::string getTypeVTK(unsigned i) const;
    
    virtual std::string getNameVTK(unsigned i) const;
    
    virtual std::vector<Mdouble> getFieldVTK(unsigned i) const;
    
    virtual void actionsAfterTimeStep()
    {};

    virtual bool isSphericalParticle() const {return false;}
    
    //const HGridCell& getHGridCell() const;
    const HGridCell& getHGridCell() const
    { return hGridCell; }

    virtual void computeMass(const ParticleSpecies& s);

protected:
    
    Mdouble radius_; ///Particle radius_
    Mdouble invMass_; ///Inverse Particle mass (for computation optimization)
    MatrixSymmetric3D invInertia_; ///Inverse Particle inverse inertia (for computation optimization)


private:
    
    /*!
     * \brief Pointer to the particle's ParticleHandler container
     */
    ParticleHandler* handler_;

#ifdef CONTACT_LIST_HGRID
    /*!
     * \brief 
     */
    PossibleContact* firstPossibleContact;
#endif
    
    /// All hGrid-information: the cell (x,y,z,level), and the previous and next particle in this cell compared to this
    /// particle
    HGridCell hGridCell;
    
    BaseParticle* hGridNextObject_; ///Pointer to the next Particle in the same HGrid cell
    BaseParticle* hGridPrevObject_; ///Pointer to the previous Particle in the same HGrid cell
    
    ///Particle attributes
    BaseParticle* periodicFromParticle_; ///Pointer to originating Particle
    
    //MPI particle attributes
    bool isMPIParticle_; /// returns true if the particle acts as an MPI particle instead of a real particle
    bool isInMPIDomain_; /// returns true if it flagged as being in MPI domain
    unsigned communicationComplexity_; //Determins if the particle in the corner of a domain, rib, side or not at all
    
    //Periodic boundary particle attributes
    bool isInPeriodicDomain_; /// bool that indicates if a particle is in the periodic domain of any boundary
    bool isPeriodicGhostParticle_; /// Indicates if the particle is a ghost particle of a periodic particle
    std::vector<int> previousPeriodicComplexity_; /// Indicates the periodic complexity at previous time step
    std::vector<int> periodicComplexity_;  /// Indicates the periodic complexity at current time step. Used to update periodic status
    bool isMaserParticle_; ///Indicates if this particle belongs to the maser boundary or is released into the wide open world.
    
    Vec3D displacement_; ///Displacement (only used in StatisticsVector, StatisticsPoint)
    Vec3D previousPosition_; /// Particle's position at previous time step
    
    /*!
     * Since ParticleSpecies is allowed to set the mass of a BaseParticle, it is
     * a friend of this class.
     */
    friend void ParticleSpecies::computeMass(BaseParticle*) const;
    
    Mdouble info_; // by default, the species ID (if a species is set)
};

#endif
