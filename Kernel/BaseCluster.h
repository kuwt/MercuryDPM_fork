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

#ifndef BaseCluster_H
#define BaseCluster_H

#endif //BaseCluster_H

#include <Mercury3D.h>
#include "Species/LinearPlasticViscoelasticFrictionSpecies.h"


#ifdef MERCURY_USE_MPI
#include <mpi.h>
#include <MpiDataClass.h>
#include <MpiContainer.h>
#endif

/*!
 * \class ClusterDPM
 * \brief An object of this class is inside FixedClusterInsertionBoundary and RandomClusterInsertionBoundary.
 * \details Inherits from Mercury3D: this is the core of the cluster generation process. An object of this class is
 *          called inside BaseClusterInsertionBoundary::checkBoundaryBeforeTimeStep and
 *          FixedClusterInsertionBoundary::checkBoundaryBeforeTimeStep.
 */
#ifndef BaseCluster_h
#define BaseCluster_h
class BaseCluster : public Mercury3D
{
public:

/*
 * ----------------------------------------------------------------------
 *               FUNCTIONS: setters and getters
 * ----------------------------------------------------------------------
 */

    /*!
     * \brief Default constructor.
     */
    BaseCluster();

    /*!
     * \brief Default destructor.
     */
    ~BaseCluster() final;

    /*!
     * \brief This returns the value of position_, which is the position in which the cluster will be inserted.
     */
    Vec3D getPosition() const;

    /*!
     * \brief This sets the value of position_, which is the position in which the cluster will be inserted.
     */
    void setPosition(Vec3D p);

    /*!
     * \brief This returns the value of the ratio between collision time and time step.
     */
    Mdouble getCollisionTimeOverTimeStep() const;

    /*!
     * \brief This sets the collisionTimeOverTimeStep number (which is the ratio between collision time and time step).
     */
    void setCollisionTimeOverTimeStep(Mdouble cTOTS);

    /*!
     * \brief This returns the value of particles' radius if there's no dispersity in size. In case of dispersity != 1,
     *          this is the radius from which all radii are computed (as a consequence in this case it is also the
     *          pseudo-averaged radius).
     */
    Mdouble getRadiusParticle() const;

    /*!
     * \brief This sets the value of particles' radius if there's no dispersity in size.
     */
    void setRadiusParticle(Mdouble rP);

    /*!
     * \brief This returns the value of particles' dispersity in size.
     */
    Mdouble getSizeDispersityParticle() const;

    /*!
     * \brief This sets the value of particles' dispersity in size.
     */
    void setSizeDispersityParticle(Mdouble sDP);

    /*!
     * \brief This returns the value of the number of particles in the cluster.
     */
    int getNumberOfParticles() const;

    /*!
     * \brief This sets the value of the number of particles in the cluster.
     */
    void setNumberOfParticles(int nP);

    /*!
     * \brief This sets the desired value of the cluster radius (there is no getter of this value, but there is a getter
 *            of the actual mean cluster radius obtained, getMeanClusterRadius)
     */
    void setRadiusCluster(Mdouble rCR);

    /*!
     * \brief This gets the final value obtained for the mass fraction;
     */
    Mdouble getFinalMassFraction();

    /*!
     * \brief This returns the value of the cluster ID.
     */
    unsigned int getClusterId() const;

    /*!
     * \brief This sets the value of the cluster ID.
     */
    void setClusterId(unsigned int iC);

    /*!
     * \brief This returns the value of the velocity damping modulus.
     */
    Mdouble getVelocityDampingModulus() const;

    /*!
     * \brief This sets the value of the velocity damping modulus.
     */
    void setVelocityDampingModulus(Mdouble vDM);

    /*!
     * \brief This returns the value of the number of particles used to compute internal structure
     */
    int getNumberOfInternalStructurePoints() const;

    /*!
     * \brief This sets the value of the number of particles used to compute the internal structure.
     */
    void setNumberOfInternalStructurePoints(int gL);

    /*!
     * \brief This returns the value of the value of the energy ratio threshold under which the process can be
     *          considered static, and so over.
     */
    Mdouble getEnergyRatioTolerance() const;

    /*!
     * \brief This sets the value of the value of the energy ratio threshold under which the process can be
     *          considered static, and so over.
     */
    void setEnergyRatioTolerance(Mdouble eRT);

    /*!
     * \brief This returns the species of the particle.
     */
    LinearPlasticViscoelasticFrictionSpecies *getParticleSpecies() const;

    /*!
     * \brief This sets the species of the particle.
     */
    void setParticleSpecies(LinearPlasticViscoelasticFrictionSpecies *particleSpecies);

    /*!
     * \brief This gets the value of velocity after creation.
     */

    Vec3D getVelocity();
    
    /*!
     * \brief This sets the value of velocity after creation.
     */

    void setVelocity(Vec3D v);

    /*!
     * \brief This returns the bool variable that defines whether the cluster data output (which is NOT the mercury
     *          data output) is written or not.
     */
    bool isCdatOutputOn() const;

    /*!
     * \brief This sets the bool variable that defines whether the cluster data output will be written or not.
     */
    void doCdatOutput(bool iCOO);

    /*!
     * \brief This returns the bool variable that defines whether the cluster overlap output is written or not.
     */
    bool isOverlOutputOn() const;

    /*!
     * \brief This sets the bool variable that defines whether the cluster overlap output will be written or not.
     */
    void doOverlOutput(bool iOOO);

    /*!
     * \brief This returns the bool variable that defines whether the cluster adjacency matrix output is written or not.
     */
    bool isAmatOutputOn() const;

    /*!
     * \brief This sets the bool variable that defines whether the cluster adjacency matrix output will be written
     *          or not.
     */
    void doAmatOutput(bool iAOO);

    /*!
     * \brief This returns the bool variable that defines whether the cluster internal structure output is written or not.
     */
    bool isIntStrucOutputOn() const;

    /*!
     * \brief This sets the bool variable that defines whether the cluster internal structure output will be written
     *          or not.
     */
    void doIntStrucOutput(bool iISOO);

    /*!
     * \brief This returns the bool variable that defines whether the cluster vtk output is written or not.
     */
    bool isVtkOutputOn() const;

    /*!
     * \brief This sets the bool variable that defines whether the cluster vtk output will be written or not.
     */
    void doVtkOutput(bool iVOO);

    /*!
     * \brief This returns the bool variable that defines whether the cluster restart output is written or not.
     */
    bool isRestartOutputOn() const;

    /*!
     * \brief This sets the bool variable that defines whether the cluster restart output will be written or not.
     */
    void doRestartOutput(bool isRestartOutputOn);

    /*!
     * \brief This returns the bool variable that defines whether the cluster fStat output is written or not.
     */
    bool isFStatOutputOn() const;

    /*!
     * \brief This sets the bool variable that defines whether the cluster fStat output will be written or not.
     */
    void doFStatOutput(bool isfStatOutputOn);

    /*!
     * \brief This returns the bool variable that defines whether the cluster ene output is written or not.
     */
    bool isEneOutputOn() const;

    /*!
     * \brief This sets the bool variable that defines whether the cluster ene output will be written or not.
     */
    void doEneOutput(bool isEneOutputOn);

    /*!
     * \brief this returns meanClusterRadius (radius of an ideal perfectly spherical cluster, there's no setter).
     */
    Mdouble getMeanClusterRadius();

    /*!
     * \brief this returns the average overlap.
     */
    Mdouble getAverageOverlap();

/*
 * ----------------------------------------------------------------------
 *               FUNCTIONS: overridden mercury3D functions
 * ----------------------------------------------------------------------
 */


    /*!
     * \brief Overrides DPMBase setupInitialConditions(): in this initial conditions for the problem are set.
     */
    void setupInitialConditions() override;

    /*!
     * \brief Overrides DPMBase actionsAfterTimeStep(): in this compression and decompression are computed,
     *          depending on the variable stage_.
     */
    void actionsAfterTimeStep() override;

    /*!
     * \brief Overrides DPMBase actionsAfterSolve(): in this cluster data file and cluster overlap file are closed and
     *          a few final actions are executed. Details in BaseCluster.cc.
     */
    void actionsAfterSolve() override;

    /*!
     * \brief Overrides DPMBase write(): in this all variables needed by the program for restarting are written.
     */
    void write(std::ostream& os, bool writeAllParticles ) const override;

    /*!
     * \brief Overrides DPMBase read(): in this all variables needed by the program for restarting are read.
     */
    void read(std::istream& is, ReadOptions opt = ReadOptions::ReadAll) override;

    /*!
     * \brief Overrides DPMBase actionsOnRestart(): in this all variables needed by the program for restarting
     *          are initialized.
     */
    void actionsOnRestart() override;

    /*!
     * \brief Overrides DPMBase printTime(): this way variables of interest are shown.
     */
    void printTime() const override;
    
private:
    /*
     * ----------------------------------------------------------------------
     *           FUNCTIONS: functions inside setupInitialConditions
     * ----------------------------------------------------------------------
     */

    /*!
     * \brief Sets all radii according to particleRadius and sizeDispersityParticle.
     */
    void setRadii();

    /*!
     * \brief Sets species of particles.
     */
    void setSpecies();

    /*!
     * \brief Sets domain limits.
     */
    void setDomainLimits();

    /*!
     * \brief Calculates the time step.
     */
    void calculateTimeStep();

    /*!
     * \brief Inserts particles inside the domain.
     */
    void insertParticles();

    /*!
     * \brief Creates the cluster data output file.
     */
    void makeCdatFile();

    /*!
     * \brief Creates the cluster overlap output file.
     */
    void makeOverlFile();

    /*!
     * \brief This function tries to insert the n-th particle (returns true if it manage to do that).
     *          It is inside insertParticles().
     */
    bool particleInsertionSuccessful(int n);



    /*
     * ----------------------------------------------------------------------
     *        FUNCTIONS: functions inside actionsAfterTimeStep
     * ----------------------------------------------------------------------
     */

    /*!
     * \brief This functions computes some important cluster information needed by the program.
     */
    void makeDataAnalysis();

    /*!
     * \brief This writes on the cluster data output file.
     */
    void writeToCdatFile();

    /*!
     * \brief This writes on the cluster overlap output file.
     */
    void writeToOverlFile();

    /*!
     * \brief This applies force on each particle.
     */
    void applyCentralForce();

    /*!
     * \brief This linearly increases the value of forceModulus (stage = 1).
     */
    void increaseForce();

    /*!
     * \brief This damps values of each particle velocity (stage = 1, stage = 2, stage = 3).
     */
    void dampVelocities();

    /*!
     * \brief This linearly decreases values of forceModulus (stage = 2).
     */
    void decreaseForce();

    /*!
     * \brief This damps values of forceModulus (stage = 3).
     */
    void dampForce();

    /*!
     * \brief This calculates the adjacency matrix of the cluster.
     */
    void createAdjacencyMatrix();

    /*!
     * \brief This creates the adjacency matrix file.
     */
    void makeAmatFile();

    /*!
     * \brief This writes on the adjacency matrix file.
     */
    void writeAmatFile();

    /*!
     * \brief This computes the internal structure of the cluster.
     */
    void computeInternalStructure();

    /*!
     * \brief This creates the gnuplot file needed for printing force vs overlaps values.
     */
    void makeGnuplotFile();

    /*!
     * \brief This creates the file needed for writing down datas from computeInternalStructure().
     */
    void makeIntenalStructureFile();



    /*
     * ----------------------------------------------------------------------
     *       VARIABLES: accessed by the user with setters and getters
     * ----------------------------------------------------------------------
     */

    //POSITION
    //\brief Position where the cluster is inserted after creation.
    Vec3D position_;

    // TIME
    //\brief Ratio between collision time and time step: should be at least 50.
    Mdouble collisionTimeOverTimeStep_;
    // \brief Energy ratio threshold under wich the simulation can be considered static and so can be stopped.
    Mdouble energyRatioTolerance_;

    // PARTiCleS
    //\brief Radius basing on which all radii will be computed.
    Mdouble radiusParticle_;
    // \brief Bool that saves whether or not the user has set the radius of the single particle composing the cluster.
    bool setRadiusParticle_ = false;
    //\brief Size dispersity of particles: must be between >= than 1.
    Mdouble sizeDispersityParticle_;
    //\brief Total number of particles.
    int nParticles_;
    // \brief Bool that saves whether or not the user has set the number of particles
    bool setNumberOfParticles_ = false;

    // CLUSTER
    //\brief Total number of particles.
    unsigned int idCluster_;
    //\brief Desired radius of the cluster.
    Mdouble radiusCluster_;
    // \brief Bool that saves whether or not the user has set the cluster radius
    bool setRadiusCluster_ = false;
    //\brief Velocity of the cluster after creation
    Vec3D clusterVelocity_;
    //\brief mean cluster radius after creation.
    Mdouble meanClusterRadius_;

    // CENTRAL FORCE
    //\brief Value of damping modulus for velocity.
    Mdouble velocityDampingModulus_;

    //  DATA ANALYSIS
    // \brief Number of points used for computing internal structure.
    int nInternalStructurePoints_;

    // SPECIES
    //\brief particle species.
    LinearPlasticViscoelasticFrictionSpecies* particleSpecies_;


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
    //\brief bool used to define whether or not vtk output must be created.
    bool isRestartOutputOn_;
    //\brief bool used to define whether or not vtk output must be created.
    bool isFStatOutputOn_;
    //\brief bool used to define whether or not vtk output must be created.
    bool isEneOutputOn_;


    /*
     * ----------------------------------------------------------------------
     *                VARIABLES: never accessed by the user
     * ----------------------------------------------------------------------
     */


    // PARTICLES
    //\brief vector in which all radii will be stored after calculation.
    std::vector<Mdouble> radii_;
    //smallest radius
    Mdouble smallestRadius_;
    //\brief mass of the particle which has radius radiusParticle.
    Mdouble massParticle_;
    //\brief total volume of all particles.
    Mdouble totalParticleVolume_;

    // GEOMETRY
    //\brief size of the cubic domain.
    Mdouble boxSize_;
    //\brief center of mass.
    Vec3D centerOfMass_;

    // CONTACT RELATED
    //\brief adjacency matrix.
    std::vector< std::vector<int> > adjacencyMatrix_;
    //\brief mean coordination number.
    Mdouble meanCoordinationNumber_;
    //\brief maximum relative overlap.
    Mdouble maxRelativeOverlap_;
    //\brief mean relative overlap.
    Mdouble meanRelativeOverlap_;
    //\brief minimum relative overlap.
    Mdouble minRelativeOverlap_;
    //\brief number of total intra-cluster bonds.
    int nIntraClusterBonds_;

    // OUTPUT
    //\brief cluster data file.
    std::ofstream cdatFile_;
    //\brief cluster overlap file.
    std::ofstream overlFile_;
    //\brief gnuplot file.
    std::ofstream gnuplotFile_;
    //\brief adjacency matrix file.
    std::ofstream amatFile_;
    //\brief internal structure file.
    std::ofstream intStructFile_;
    //\brief output time of files and print time.
    Mdouble fileOutputTimeInterval_;

    // DATA ANALYSIS
    //\brief radius with which solid fraction is computed
    Mdouble radiusForSolidFraction_;
    //\brief solid fraction computed with the total particle volume and radiusForSolidFraction_.
    Mdouble solidFraction_;
    //\brief solid fraction computed with internal structure analysis.
    Mdouble solidFractionIntStruct_;

    // TIME
    // stage of the simulation: 1 compression, 2 decompression, 3 relaxation, 4 simulation ended.
    int stage_;
    //\brief time flag used for the stages duration.
    Mdouble t0_;
    //\brief final time.
    Mdouble clusterTimeMax_;

    // CENTRAL FORCE
    //\brief maximum force modulus applied on particles (this value is then multiplied by distance from force center).
    Mdouble maximumForceModulus_;
    //\brief force modulus applied on particles at a certain simulation time.
    Mdouble forceModulus_;
    //\brief time interval on which force is tuned (increased or decreased).
    Mdouble forceTuningInterval_;
    //\brief time interval on which velocity is tuned (increased or decreased).
    Mdouble velocityDampingInterval_;
    //\brief time duration of force tuning (i.e. duration of compression and decompression stages).
    Mdouble forceTuningDuration_;
    //\brief maximum possible time duration of dissipation (i.e. duration of dissipation if energy ratio
    //          tollerance not reached).
    Mdouble dissipationDuration_;
    //\brief Value of damping modulus for forceFactor.
    Mdouble forceDampingModulus_;

};
#endif