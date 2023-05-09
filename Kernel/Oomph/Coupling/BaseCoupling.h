//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://MercuryDPm.org/Team>.
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

#ifndef BASE_COUPLING_H
#define BASE_COUPLING_H

#include "Math/Vector.h"
#include "Walls/TriangleWall.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "CG/Functions/Gauss.h"
#include "CG/Functions/Lucy.h"
#include "CG/Functions/Heaviside.h"
#include "Interactions/BaseInteraction.h"
#include "Logger.h"

/**
 * Common functionality for both surface- and volume-coupled problems:
 *  - defines a cg function
 *  - getParticlesInCell to find particles in a particular region
 *  - modifies the ene files to add coupled mass and energy
 *  - solveOomph and solveMercury functions to advance the two solutions in time
 *  \author Hongyang Cheng <chyalexcheng@gmail.com>
 */
template<class M, class O>
class BaseCoupling : public M, public O
{
//protected:
//    // A reference to the Mercury problem
//    Mercury3D& m = this;
//    // A reference to the Oomph problem
//    SolidProblem<typename OomphProblem::ELEMENT_TYPE>& o = this;

public:
    
    BaseCoupling() = default;
    
    void setName(std::string name) {
        M::setName(name);
        O::setName(name);
    }

    std::string getName() const {
        return M::getName();
    }

    void removeOldFiles() const {
        M::removeOldFiles();
        O::removeOldFiles();
    }
    
    /**
     * override writeEneTimeStep because mass and elastic energy are computed different for the coupling; also adds momentum and angular momentum
     */
    void writeEneTimeStep(std::ostream& os) const override
    {
        if (M::eneFile.getCounter() == 1 || M::eneFile.getFileType() == FileType::MULTIPLE_FILES ||
                M::eneFile.getFileType() == FileType::MULTIPLE_FILES_PADDED)
        {
            writeEneHeader(os);
        }

        const Mdouble m = M::particleHandler.getMass();
        const Vec3D com = M::particleHandler.getMassTimesPosition();
        //Ensure the numbers fit into a constant width column: for this we need the precision given by the operating system,
        //plus a few extra characters for characters like a minus and scientific notation.
        const static long int width = os.precision() + 6;
        os << std::setw(width) << M::getTime()
//           << " " << std::setw(width) << getCoupledMass()
           << " " << std::setw(width) << M::particleHandler.getMomentum().getX()
           << " " << std::setw(width) << M::particleHandler.getMomentum().getY()
           << " " << std::setw(width) << M::particleHandler.getMomentum().getZ()
           << " " << std::setw(width) << M::particleHandler.getAngularMomentum().getX()
           << " " << std::setw(width) << M::particleHandler.getAngularMomentum().getY()
           << " " << std::setw(width) << M::particleHandler.getAngularMomentum().getZ()
           << " " << std::setw(width) << -Vec3D::dot(M::getGravity(), com)
           << " " << std::setw(width) << M::particleHandler.getKineticEnergy()
           << " " << std::setw(width) << M::particleHandler.getRotationalEnergy()
           //<< " " << std::setw(width) << M::getElasticEnergyCoupled()
           // we need to write x, y and z coordinates separately, otherwise the width of the columns is incorrect
           << " " << std::setw(width)
           << ( m == 0 ? constants::NaN : com.X / m ) //set to nan because 0/0 implementation in gcc and clang differs
           << " " << std::setw(width) << ( m == 0 ? constants::NaN : com.Y / m )
           << " " << std::setw(width) << ( m == 0 ? constants::NaN : com.Z / m )
           << std::endl;
    }
    
//    /**
//     * elastic energy is multiplied by the coupling weight (used in writeEneTimeStep)
//     * \todo this only works for one species; rather, we should override elastic energy
//     */
//    Mdouble getElasticEnergyCoupled() const
//    {
//        // get elastic energy
//        Mdouble elasticEnergy = 0.0;
//        for (const BaseInteraction* c : MercuryProblem::interactionHandler)
//        {
//            Mdouble energyNormal = c->getNormalElasticEnergy();
//            Mdouble energyTangential = 0;
//            elasticEnergy += c->getCouplingWeight() * ( energyNormal + energyTangential );
//        }
//        return elasticEnergy;
//    }
    
//    /**
//     * If coupled, mass is multiplied by coupling weight at the center of mass (used in writeEneTimeStep)
//     */
//    Mdouble getCoupledMass() const
//    {
//        Mdouble m = 0;
//        for (BaseParticle* p : MercuryProblem::particleHandler)
//        {
//            if (!( p->isFixed() || p->isMPIParticle() || p->isPeriodicGhostParticle() ))
//            {
//                m += p->getMass() / p->getInvCouplingWeight();
//            }
//        }
//        return m;
//    }
    
    /**
     *  override writeEneHeader in DPMBase class for the coupling
     */
    void writeEneHeader(std::ostream& os) const override
    {
        //only write if we don't restart
        if (M::getAppend()) {
            return;
        }
        
        long width = os.precision() + 6;
        os << std::setw(width)
           << "time " << std::setw(width)
           << "mass " << std::setw(width)
           << "momentum_X " << std::setw(width)
           << "momentum_Y " << std::setw(width)
           << "momentum_Z " << std::setw(width)
           << "angMomentum_X " << std::setw(width)
           << "angMomentum_Y " << std::setw(width)
           << "angMomentum_Z " << std::setw(width)
           << "gravitEnergy " << std::setw(width) //gravitational potential energy
           << "traKineticEnergy " << std::setw(width) //translational kinetic energy
           << "rotKineticEnergy " << std::setw(width) //rotational kE
           << "elasticEnergy " << std::setw(width)
           << "centerOfMassX " << std::setw(width)
           << "centerOfMassY " << std::setw(width)
           << "centerOfMassZ\n";
    }
    
    /**
     * used in OomphMercuryCoupling::computeOneTimeStepForV/SCoupling to do a oomph timestep (V/SCoupling)
     */
    void solveOomph()
    {
        O::actionsBeforeOomphTimeStep();
        this->unsteady_newton_solve(this->time_pt()->dt());
    }
    
    /**
     * solve a given number of time steps nt in Mercury
     */
    void solveMercury(unsigned long nt)
    {
        for (int n = 0; n < nt; ++n)
        {
            M::computeOneTimeStep();
        }
    }
    
    /**
     * find all particles in a particular region
     * \param min
     * \param max
     * \param[out] pList used to return the list of particles
     */
    void getParticlesInCell(const Vec3D& min, const Vec3D& max, std::vector<BaseParticle*>& pList)
    {
        HGrid* const hGrid = M::getHGrid();
        
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
            if (( occupiedLevelsMask & 1 ) == 0)
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
                                pList.push_back(p);
                            }
                            p = p->getHGridNextObject();
                        }
                    }
                }
            }
        } //end for level
    }
    
    inline void setCGWidth(const double& width)
    {
        if (width == 0)
        {
            CGMapping_ = false;
        }
        else
        {
            CGMapping_ = true;
            // note CG function needs to be nondimensionalized
            // fixme CG width is not set with respect to particle radius. Should we do that?
            cgFunction_.setWidth(width);
        }
    }

    /**
     * get the dimensionalized CG width to be used in MercuryProblem
     */
    inline double getCGWidth()
    {
        return cgFunction_.getWidth();
    }
    
    inline bool useCGMapping()
    { return CGMapping_; }
    
    inline CGFunctions::LucyXYZ getCGFunction()
    { return cgFunction_; }

private:
    // flag to construct mapping with FEM basis functions (default) or DPM coarse graining
    bool CGMapping_ = false;
    // coarse-graining functions and coordinates
    CGFunctions::LucyXYZ cgFunction_;
};

#endif //BASE_COUPLING_H
