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

#ifndef SPHERICALINDENTER_H
#define	SPHERICALINDENTER_H

#include "Mercury3D.h"
//#include <assert.h>
#include "Logger.h"

/** Puts a spherical indenter into Mercury3D */
class SphericalIndenter : public Mercury3D
{
public:

    /** initialises indenter properties */
    SphericalIndenter(Mdouble indenterDiameter, Mdouble indentationVelocity, Mdouble indentationForce)
    : indenterDiameter_(indenterDiameter), indentationVelocity_(indentationVelocity),
    indentationForce_(indentationForce)
    {
        std::cout << "Creating indenter with"
            << " diameter " << indenterDiameter_
            << " force " << indentationForce_
            << " velocity " << indentationVelocity_
            << std::endl;
        logger.assert_always(indenterDiameter_ >= 0.0,"");
        logger.assert_always(indentationVelocity_ >= 0.0,"");
        logger.assert_always(indentationForce_ >= 0.0,"");
    }

    /** first create particles and species before calling this */
    void setupInitialConditions() override
    {
        std::cout << "Setting up indenter\n";
        
        indenter_.setHandler(&particleHandler);
        indenter_.setRadius(0.5 * indenterDiameter_);
        indenter_.setSpecies(speciesHandler.getLastObject());
        indenter_.fixParticle();
        indenter_.setIndex(particleHandler.getNumberOfObjects());
        indenter_.setId(particleHandler.getNumberOfObjects());
        setIndenterVelocity(-10.0*indentationVelocity_);
        logger(INFO,"Indentation Velocity %",getIndenterVelocity());

        setIndenterHeight(getBedHeight());
        
        std::cout << indenter_ << "\n";
        std::cout << "H" << getIndenterHeight() << "\n";
        std::cout << "H" << indenter_.getIndex() << "\n";
        std::cout << "H" << indenter_.getId() << "\n";
    }

    Mdouble getIndentationForce() const
    {
        return indentationForce_;
    }
    
    void setIndentationForce(Mdouble indentationForce)
    {
        indentationForce_=indentationForce;
    }

    Mdouble getIndenterHeight() const
    {
        return indenter_.getPosition().Z - indenter_.getRadius();
    }
    
    void setIndenterHeight(Mdouble height)
    {
        indenter_.setPosition(Vec3D(0.5 * (getXMax() + getXMin()),
                                    0.5 * (getYMax() + getYMin()),
                                    height+indenter_.getRadius()));
    }

    Mdouble getForceOnIndenter() const
    {
        return indenter_.getForce().Z;
    }

    void setIndenterVelocity(Mdouble indentationVelocity_)
    {
        indenter_.setVelocity({0, 0, indentationVelocity_});
    }

    Mdouble getIndenterVelocity() const
    {
        return indenter_.getVelocity().Z;
    }
    
    /** add indenter forces */
    void computeExternalForces(BaseParticle* p) override
    {
        DPMBase::computeExternalForces(p);
        if (!p->isFixed())
            computeInternalForce(p, &indenter_);
    }

    /** add indenter forces */
    void actionsBeforeTimeStep() override
    {
        indenter_.resetForceTorque(getNumberOfOMPThreads());
    }

    /** add indenter to xballs */
    void outputXBallsData(std::ostream& os) const override
    {
        // adds one line to the particle data
        os << particleHandler.getNumberOfObjects() + 1 << " " << getTime() << " "
            << getXMin() << " " << getYMin() << " " << getZMin() << " "
            << getXMax() << " " << getYMax() << " " << getZMax() << " " << std::endl;
        // This outputs the particle data
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            outputXBallsDataParticle(i, 14, os);

        os  << indenter_.getPosition() << " "
            << indenter_.getVelocity() << " "
            << indenter_.getRadius() << "  0 0 0 0 0 0 0\n";
    }

    void actionsBeforeTimeLoop() override 
    {
        std::cout << std::setw(11) << "time"
            << "\t" << std::setw(11)<< "displacement"
            << "\t" << std::setw(8) << "relForce"
            << "\t" << std::setw(9) << "direction"
            << "\t" << std::setw(11) << "eneRatio"
            << std::endl;
    }
    /** creates custom console output */
    void printTime() const override
    {
        //writeEneTimeStep(std::cout);
        std::cout << std::setw(11) << getTime()
            << "\t" << std::setw(11)<< getIndenterHeight()
            << "\t" << std::setw(11) << getForceOnIndenter() / indentationForce_
            << "\t" << std::setw(9) << (getIndenterVelocity()<0?"down":(getIndenterVelocity()==0?"stop":"up"))
            << "\t" << std::setw(11) << getKineticEnergy()/getElasticEnergy()
            << std::endl;
    }

    /** creates custom ene header */
    void writeEneHeader(std::ostream& os) const override
    {
        os << "time\tdisplacement\tforce\tdirection\n";
    }

    /** creates custom ene output */
    void writeEneTimeStep(std::ostream& os) const override
    {
        os << std::setw(12) << getTime()
            << "\t" << std::setw(12)<< getIndenterHeight()
            << "\t" << std::setw(12) << getForceOnIndenter()/ indentationForce_
            << "\t" << (getIndenterVelocity()<0?"down":(getIndenterVelocity()==0?"stop":"up"))
            << std::endl;
//        os << "t " << getTime()
//            << " \t " << (getIndenterVelocity()<0?"down":"up")
//            << "\tz " << getIndenterHeight()
//            << "\trelForce " << getForceOnIndenter() / indentationForce_ << std::endl;
    }

    void actionsAfterTimeStep() override
    {
        indenter_.sumForceTorqueOMP();
        
        static Mdouble timeToRetract = 0;
        indenter_.move(indenter_.getVelocity() * getTimeStep());

        if (getIndenterVelocity() < 0)
        { //moving down
            if (getForceOnIndenter() > 0 && getIndenterVelocity()!=-indentationVelocity_)
            {
                setIndenterVelocity(-indentationVelocity_);
                logger(INFO,"Lowering indentation velocity %",getIndenterVelocity());
            }
            if (getForceOnIndenter() >= indentationForce_)
            {
                setIndenterVelocity(0);
                timeToRetract = 1.1*getTime();
                std::cout << "stopping indenter" << std::endl;
            }
        }
        else if (getIndenterVelocity() == 0) 
        { //stopping
            if (getTime()>timeToRetract) {
                setIndenterVelocity(indentationVelocity_);
                std::cout << "retracting indenter" << std::endl;
                measuredIndentationForce = getForceOnIndenter();
                measuredIndentation = getIndenterHeight();
            }        
        } else 
        { //retracting
            if (measuredForceGradient==0 && getForceOnIndenter() < 0.5*measuredIndentationForce)
            {
                measuredForceGradient = (0.5*measuredIndentationForce)/(getIndenterHeight()-measuredIndentation);
                std::cout << "measured force gradient = " << measuredForceGradient << std::endl;
            } else if (getForceOnIndenter() > 0) {
                //std::cout << "stopping indenter" << std::endl;
                setTimeMax(getTime()*1.2);
            } else if (measuredElasticDisplacement==0){
                measuredElasticDisplacement = getIndenterHeight()-measuredIndentation;
                std::cout << "measured elastic displacement = " << measuredElasticDisplacement << std::endl;
                std::cout << "measured Elasticity = " << 0.5*measuredForceGradient/sqrt(0.5*indenterDiameter_*measuredElasticDisplacement) << std::endl;
                setIndenterVelocity(10.0*indentationVelocity_);
                logger(INFO,"Increasing indentation velocity %",getIndenterVelocity());
            }
        }
    }

    /** top of the particles (can be used as initial height of the indenter)) */
    Mdouble getBedHeight()
    {
        Mdouble bedHeight = getZMin();
        for (auto p : particleHandler)
            if (!p->isFixed())
                bedHeight = std::max(bedHeight, p->getPosition().Z + p->getRadius());
        std::cout << "bed height" << bedHeight << " N=" << particleHandler.getNumberOfObjects() << " \n";
        //return 1;
        return bedHeight;
    }

private:
    SphericalParticle indenter_;
    Mdouble indenterDiameter_;
    Mdouble indentationVelocity_;
    Mdouble indentationForce_;
    
    Mdouble measuredIndentationForce = 0;
    Mdouble measuredIndentation = 0;
    Mdouble measuredForceGradient = 0;
    Mdouble measuredElasticDisplacement = 0;
    
};

#endif	/* SPHERICALINDENTER_H */

