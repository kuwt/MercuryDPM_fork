//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include "Mercury2D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/SphericalParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "CG/TimeAveragedCG.h"
#include "CG/CG.h"

enum BottomType
{
    FLAT, ORDERED, DISORDERED
};

class ChutePeriodic2D : public Mercury2D
{
public:

    /* constructor; default parameters */

    ChutePeriodic2D() //assuming domain has a lower corner (0,0,0); YMin may extend if rough bottom
    {
        setName("ChutePeriodic2D");

        //Set and add the particle-species.
        LinearViscoelasticFrictionSpecies species;
        species.setDensity(4 / constants::pi);
        species.setCollisionTimeAndRestitutionCoefficient(5e-3, 0.88, 1);
        species.setSlidingDissipation(species.getDissipation());
        species.setSlidingStiffness(2. / 7. * species.getStiffness());
        species.setSlidingFrictionCoefficient(0.5);
        speciesHandler.copyAndAddObject(species);
        
        // Chute geometry
        setChuteAngleAndMagnitudeOfGravity(22,9.81);
        setChuteLength(50);
        setInflowHeight(15);

        // Particle size
        setMinInflowParticleRadius(0.4);
        setMaxInflowParticleRadius(0.6);
        setBottomType(DISORDERED);
        setFixedParticleRadiusRatio(1.2);
        setFixedParticleSpacing(0.5);

        //Set time control parameters
        setTimeMax(20);
        setTimeStep(1.0e-4); //(1/50th of the collision time)
        setSaveCount(10000);
        setIsFlowing(true);

        //Set output files
        eneFile.setFileType(FileType::ONE_FILE);
        dataFile.setFileType(FileType::NO_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::NO_FILE);

        //Set CG properties
        CG<CGCoordinates::Y> cg0;
        cg0.statFile.setSaveCount(dataFile.getSaveCount());
        cg0.setY(getYMin(),getYMax());
        cg0.setHY(getInflowParticleRadius());
        cg0.setZ(0,1);
        cgHandler.copyAndAddObject(cg0);
        TimeAveragedCG<CGCoordinates::Y> cg1;
        cg1.setY(getYMin(),getYMax());
        cg1.setHY(getInflowParticleRadius());
        cg1.setZ(0,1);
        cg1.setTimeMin(5);
        cgHandler.copyAndAddObject(cg1);
    }

    /* some checks/updates after each time step */

    void actionsAfterTimeStep() override
    {
        DPMBase::actionsAfterTimeStep();

//        setYMax(particleHandler.getHighestPositionComponentParticle(1)->getPosition().Y);
//        cgHandler.getObject(0)->setY(getYMin(),getYMax());
//        cgHandler.getObject(1)->setY(getYMin(),getYMax());

        // check flow state every 10 time units
        if (getTime() > 100)
        {
            // update flow state and print out
            if (getKineticEnergy() / getElasticEnergy() < 1e-2)
                setIsFlowing(false);
        }
    }

    bool continueSolve() const override
    {
        DPMBase::continueSolve();
        return (isFlowing_);
    }

    /* set up initial conditions */

    void setupInitialConditions() override
    {
        DPMBase::setupInitialConditions();

        // make chute periodic in x-direction:
        makeChutePeriodic();

        // creates the bottom of the chute
        createBottom();

        // create flow particles
        createFlowParticles();
    }

    void createFlowParticles()
    {
        logger(INFO, "Creating flow particles");

        // create one particle
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));

        // place/copy to different locations
        auto Nx = (int) floor(getChuteLength() / getInflowParticleRadius() / 2);
        auto Ny = (int) floor(getInflowHeight() / getInflowParticleRadius() / 2);
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {

                p0.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
                p0.setPosition(Vec3D(getInflowParticleRadius() * (2 * i + 1),
                                     getInflowParticleRadius() * (2 * j + 1),
                                     0.0));
                p0.setVelocity(Vec3D(random.getRandomNumber(-1,1) * getInflowParticleRadius() * 0.01,
                                     random.getRandomNumber(-1,1) * getInflowParticleRadius() * 0.01,
                                     0.0));
                particleHandler.copyAndAddObject(p0);
            }
        }
    }

    void createBottom()
    {

        // create fixed particles
        if (getBottomType() == ORDERED || getBottomType() == DISORDERED)
        {
            // create one fixed particles
            SphericalParticle b0;
            b0.setSpecies(speciesHandler.getObject(0));
            b0.setRadius(getFixedParticleRadius());
            // place/copy to different locations
            auto Nb = (int) floor((getChuteLength()/getFixedParticleRadius()/2-1)/(1+getFixedParticleSpacing())+1);
            for (int i = 0; i < Nb; i++) {
                Mdouble currentX = getFixedParticleRadius() * (2 * (1 + getFixedParticleSpacing()) * i + 1);
                Mdouble currentY = -getFixedParticleRadius();
                if (getBottomType() == DISORDERED) // a random offset between [-0.5,0.5]*\eps*d_b
                    currentX += random.getRandomNumber(-1,1) * getFixedParticleRadius() * getFixedParticleSpacing();
                b0.setPosition(Vec3D(currentX,currentY,0.0));
                b0.fixParticle();
                particleHandler.copyAndAddObject(b0);
            }
        }

        // create bottom wall (Y-direction)
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
    }

    void makeChutePeriodic()
    {
        //make the chute periodic in X-direction:
        PeriodicBoundary b0;
        b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
    }

    /* setters and getters */

    void setChuteAngle(Mdouble chuteAngle)
    {
        // retrieve the magnitude of gravity
        Mdouble gravity = getGravity().getLength();
        if (gravity == 0)
        {
            logger(WARN, "[Chute::setChuteAngle()] zero gravity");
        }

        // reset the gravity vector, with the given angle
        setChuteAngleAndMagnitudeOfGravity(chuteAngle, gravity);

    }

    Mdouble getChuteAngle() const
    {
        return chuteAngle_;
    }

    void setChuteAngleAndMagnitudeOfGravity(Mdouble chuteAngle, Mdouble gravity)
    {
        chuteAngle_ = chuteAngle * constants::pi / 180.0;
        setGravity(Vec3D(sin(chuteAngle_), -cos(chuteAngle_), 0.0) * gravity);
    }

    void setChuteLength(Mdouble chuteLength)
    {
        setXMax(chuteLength);
    }

    Mdouble getChuteLength() const
    {
        return getXMax();
    }

    void setInflowHeight(Mdouble inflowHeight)
    {
        inflowHeight_ = inflowHeight;
        // flow height is Y-direction
        setYMax(1.2 * inflowHeight_);
    }

    Mdouble getInflowHeight() const
    {
        return inflowHeight_;
    }

    void setInflowParticleRadius(Mdouble inflowParticleRadius)
    {
        minInflowParticleRadius_ = inflowParticleRadius;
        maxInflowParticleRadius_ = inflowParticleRadius;
    }

    Mdouble getInflowParticleRadius() const
    {
        return 0.5 * (minInflowParticleRadius_ + maxInflowParticleRadius_);
    }

    void setMinInflowParticleRadius(Mdouble minInflowParticleRadius)
    {
        minInflowParticleRadius_ = minInflowParticleRadius;
    }

    Mdouble getMinInflowParticleRadius() const
    {
        return minInflowParticleRadius_;
    }

    void setMaxInflowParticleRadius(Mdouble maxInflowParticleRadius)
    {
        maxInflowParticleRadius_ = maxInflowParticleRadius;
    }

    Mdouble getMaxInflowParticleRadius() const
    {
        return maxInflowParticleRadius_;
    }

    void setFixedParticleRadius(Mdouble fixedParticleRadius)
    {
        if (fixedParticleRadius < 1e-6)
            fixedParticleRadius = 0.0;
        fixedParticleRadius_ = fixedParticleRadius;
        setYMin(fixedParticleRadius_ * -2.0);
    }

    Mdouble getFixedParticleRadius() const
    {
        return fixedParticleRadius_;
    }

    void setFixedParticleRadiusRatio(Mdouble fixedParticleRadiusRatio)
    {
        fixedParticleRadiusRatio_ = fixedParticleRadiusRatio;
        setFixedParticleRadius(getInflowParticleRadius()*fixedParticleRadiusRatio_);
    }

    Mdouble getFixedParticleRadiusRatio() const
    {
        return fixedParticleRadiusRatio_;
    }


    void setFixedParticleSpacing(Mdouble fixedParticleSpacing)
    {
        fixedParticleSpacing_ = fixedParticleSpacing;
    }

    Mdouble getFixedParticleSpacing() const
    {
        return fixedParticleSpacing_;
    }

    void setBottomType(BottomType bottomType)
    {
        bottomType_ = bottomType;
        if (bottomType == FLAT)
            setFixedParticleRadius(0.0);
        setYMin(getFixedParticleRadius() * -2.0);
    }

    BottomType getBottomType() const
    {
        return bottomType_;
    }

    bool getIsFlowing() const
    {
        return isFlowing_;
    }

    void setIsFlowing(bool isFlowing)
    {
        isFlowing_ = isFlowing;
    }

private:

    /*!
     * \brief chute angle in degrees
     */
    Mdouble chuteAngle_;
    /*!
     * \brief radius of the fixed particles at the bottom
     */
    Mdouble fixedParticleRadius_;
    /*!
     * \brief size ratio of fixed particles to flow particles
     */
    Mdouble fixedParticleRadiusRatio_;
    /*!
    * \brief spacing of the fixed particles at the bottom (center distance / diameter - 1)
    */
    Mdouble fixedParticleSpacing_;
    /*!
     * \brief minimal radius of inflowing particles
     */
    Mdouble minInflowParticleRadius_;
    /*!
     * \brief maximal radius of inflowing particles
     */
    Mdouble maxInflowParticleRadius_;
    /*!
     * \brief Height of inflow
     */
    Mdouble inflowHeight_;
    /*!
     * \brief Determines the type of rough bottom created (if any). See also the enum
     * RoughBottomType at the beginning of this header file.
     */
    BottomType bottomType_;
    /*!
     * \brief Flowing or not
     */
    bool isFlowing_;
};

