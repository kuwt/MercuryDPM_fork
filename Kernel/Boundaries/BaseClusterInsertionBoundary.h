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

#ifndef BOUNDARIES_BASECLUSTERINSERTIONBOUNDARY_H
#define BOUNDARIES_BASECLUSTERINSERTIONBOUNDARY_H



#include "InsertionBoundary.h"
#include "BaseCluster.h"
#include "Math/Vector.h"
#include "Math/PSD.h"

class BaseParticle;

class RNG;

/*!
 * \class ClusterInsertionBoundary
 * \brief It's an insertion boundary which has cuboidal shape and inserts clusters.
 *          Two classes (RandomClusterInsertionBoundary and FixedClusterInsertionBoundary) derive from this.
 */

class BaseClusterInsertionBoundary : public InsertionBoundary
{
public:
    
    /*!
     * \brief Constructor; sets everything to 0.
     */
    BaseClusterInsertionBoundary();

    /*!
     * \brief Copy constructor with deep copy.
     */
    BaseClusterInsertionBoundary(const BaseClusterInsertionBoundary& other);

    /*!
     * \brief Destructor: default destructor.
     */
    ~BaseClusterInsertionBoundary() override;

    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    BaseClusterInsertionBoundary* copy() const override;

    //!\brief this turns off the randomise(): created for UnitTests.
    void setRandomised(bool randomised);

    //!\brief this returns a bool which indicates if the process is randomised (true) or not (false).
    bool getRandomised();

    //!\brief this returns the number of cluster inserted.
    unsigned int getNumberOfClusterInserted();

    //!\brief this sets the radius of the micro particle composing the cluster.
    void setRadiusMicroParticle(Mdouble rMP);

    /*!
    * \brief Sets the range of cluster radius that may be generated.
    */
    void setRadiusRange(Mdouble radMin, Mdouble radMax);
    
    /*!
     * \brief Sets the geometry (position and velocity distribution) of the
     * ClusterInsertionBoundary
     */
    void setGeometry(Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax);

    /*!
     * \brief Sets the velocity range of the ClusterInsertionBoundary
     */
    void setVelocityRange(Vec3D velMin, Vec3D velMax);

    /*!
     * \brief sets additional cluster properties as:
     * - collision time over time step ratio
     * - velocity damping modulus
     * - energy ratio tolerance
     */
    void setAdditionalClusterProperties(Mdouble collisionTimeOverTimeStep, Mdouble velocityDampingModulus, Mdouble energyRatioTolerance);

    /*!
     * \brief sets cluster whether or not cluster output files will be created, for example:
     * - cluster data
     * - overlap output
     * - adjacensy matrix output
     * - internal structure output
     * - vtk output
     * - restart file output
     * - fstat file output
     * - ene file output
     * For more information and details see class BaseCluster.
     */
    void setOutputClusterProperties(bool doCdatOutput, bool doOverlOutput, bool doAmatOutput, bool doIntStrucOutput,
                                    bool doVtkOutput, bool doRestartOutput, bool doFStatOutput, bool doEneOutput);

    /*!
     * \brief Generates a random position, velocity for the cluster p
     */
    void placeParticle(BaseParticle* p, RNG& random) override;

    /*!
     * \brief Fills the boundary with clusters.
     */
    void checkBoundaryBeforeTimeStep(DPMBase* md) override;

    /*!
     * \brief reads boundary properties from istream
     */
    void read(std::istream& is) override;

    /*!
     * \brief writes boundary properties to ostream
     */
    void write(std::ostream& os) const override;


private:

    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;

    /*
    * ----------------------------------------------------
    *                    VARIABLES
    * ----------------------------------------------------
    */

    // Position in which the cluster will be created and also added in the domain.
    Vec3D position_;


protected:

    //number of cluster inserted
    unsigned int nClusterInserted_;

    //Radius of the particle composing the cluster
    Mdouble radiusParticle_;

    //Particles
    //\brief Total number of particles.
    int nParticles_;

    //species with which the cluster will be created.
    LinearPlasticViscoelasticFrictionSpecies* clusterSpecies_;

    // File
    //\brief bool used to define whether or not cluster data output must be created.
    bool isCdatOutputOn_;
    //\brief bool used to define whether or not overlap data output must be created.
    bool isOverlOutputOn_;
    //\brief bool used to define whether or not adjacency matrix output must be created.
    bool isAmatOutputOn_;
    //\brief bool used to define whether or not cluster internal structure output must be created.
    bool isIntStrucOutputOn_;
    //\brief bool used to define whether or not vtk output must be created.
    bool isVtkOutputOn_;
    //\brief bool used to define whether or not restart output must be created.
    bool isRestartOutputOn_;
    //\brief bool used to define whether or not fStat output must be created.
    bool isFStatOutputOn_;
    //\brief bool used to define whether or not eneOutput output must be created.
    bool isEneOutputOn_;


    //\brief Size dispersity of particles: must be between 0 and 1
    Mdouble sizeDispersityParticle_;
    //\brief Value of damping modulus for velocity.
    Mdouble velocityDampingModulus_;
    // \brief Number of points used for creating internal structure's file.
    int nInternalStructurePoints_;
    // \brief Energy ratio threshold under wich the simulation can be considered static.
    Mdouble energyRatioTolerance_;

    //\brief Ratio between collision time and time step: should be at least 50.
    Mdouble collisionTimeOverTimeStep_;

    /*
     * \brief Minimal and maximal positions defining the boundary's boundaries,
     * and minimum and maximum velocity of the particles to be inserted.
     */
    Vec3D posMin_, posMax_, velMin_, velMax_;

    // Bool defining if the user has set the radius of a single particle composing the cluster OR the number of
    // particles inside the cluster.
    bool setRadiusParticleAndNotNumberOfParticles_;

    // Vectors defining cluster radii and position for class FixedClusterInsertionBoundary
     std::vector<Vec3D> clusterPositions_;
     std::vector<Mdouble> clusterRadii_;

    //Variable used to switch the randomise() process
    bool randomised_;

};

/*
 * write to file
 */
std::ostream& operator<<(std::ostream& os, BaseClusterInsertionBoundary::Distribution type);

/*
 * read from file
 */
std::istream& operator>>(std::istream& is, BaseClusterInsertionBoundary::Distribution& type);

#endif

