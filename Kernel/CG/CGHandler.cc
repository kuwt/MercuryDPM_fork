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
#include "CG/CGHandler.h"
#include "DPMBase.h"
#include "Walls/InfiniteWall.h"

CGHandler::CGHandler()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "CGHandler::CGHandler() finished" << std::endl;
#endif
}

/*!
 * \param[in] ch The CGHandler that has to be copied.
 */
CGHandler::CGHandler(const CGHandler& ch)
{
    setDPMBase(ch.getDPMBase());
    //copyContentsFromOtherHandler(ch);
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "CGHandler::CGHandler(const CGHandler&) finished" << std::endl;
#endif
}

/*!
 * \param[in] rhs The WallHandler on the right hand side of the assignment.
 * \details This is not a copy assignment operator! It only copies the pointer to the
 *          DPMBase and the BaseWall in objects_, it sets the other data members
 *          to 0 or nullptr.
 */
CGHandler& CGHandler::operator=(const CGHandler& rhs)
{
    if (this != &rhs)
    {
        clear();
        setDPMBase(rhs.getDPMBase());
        copyContentsFromOtherHandler(rhs);
    }
    return *this;
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "CGHandler::operator =(const CGHandler&) finished" << std::endl;
#endif
}

CGHandler::~CGHandler()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "CGHandler::~CGHandler() finished" << std::endl;
#endif
}

/*!
 * \param[in] cg A pointer to the CG opject that has to be added to the handler.
 */
void CGHandler::addObject(BaseCG* cg)
{
    //Puts the cg object in the list
    BaseHandler<BaseCG>::addObject(cg);
    //set the CGHandler pointer
    cg->setHandler(this);
}

/*!
 * \return a string containing the name.
 */
std::string CGHandler::getName() const
{
    return "CGHandler";
}

void CGHandler::readAndAddObject(std::istream& is UNUSED)
{

}

void CGHandler::write(std::ostream& os UNUSED) const
{

}

void CGHandler::initialise()
{
    //initialise all CG objects in the handler
    for (BaseCG* it : *this)
        it->initialise();
};

void CGHandler::evaluate()
{
    //evaluate all CG objects in the handler
    for (BaseCG* it : *this)
    {
        //if we are below timeMax and the next time step should be written
        if (it->getTimeMin() <= getDPMBase()->getTime() && it->getTimeMax() > getDPMBase()->getTime()
            && it->statFile.saveCurrentTimestep(getDPMBase()->getNtimeSteps()))
        {
            //logger(INFO,"evaluate %, nt=%",it->statFile.getName(), getDPMBase()->getNtimeSteps());
            it->statFile.setLastSavedTimeStep(getDPMBase()->getNtimeSteps());
            it->evaluate();
        }
    }
};

void CGHandler::finish()
{
    //enforce that data gets written
    for (BaseCG* it : *this)
    {
        it->statFile.setLastSavedTimeStep(NEVER);
    }
    //evaluate all CG objects in the handler
    //evaluate();
    //finish all CG objects in the handler
    for (BaseCG* it : *this)
        it->finish();
}

void CGHandler::restart(std::string name)
{

    // determines if the variable name contains the name of a restart file, or a problem name
    // (in which case the restart file is assumed to be name.restart)
    DPMBase* dpm = getDPMBase();
    std::string::size_type endName = name.find(".restart");
    if (endName == std::string::npos)
    {
        dpm->setName(name);
    } else
    {
        dpm->setName(name.substr(0, endName));
        dpm->restartFile.setName(name);
    }

    // set name
    cgLogger(INFO, "Evaluating files named %", dpm->getName());

    // read restart file
    if (!dpm->readRestartFile())
    {
        //what to do if restart file could not be loaded
        cgLogger(INFO, "Using default dpm setup as % does not exist", dpm->restartFile.getName());
    } else
    {
        cgLogger(INFO, "Successfully restarted from %, t=%, Np=%, Nc=%", dpm->restartFile.getFullName(), dpm->getTime(),
                 dpm->particleHandler.getSize(), dpm->interactionHandler.getNumberOfObjects());
    }
}

void CGHandler::restartAndEvaluateRestartFiles(const std::string& name)
{
    restart(name);
    evaluateRestartFiles();
}

void CGHandler::restartAndEvaluateDataFiles(const std::string& name, bool evaluateFStatFiles)
{
    restart(name);
    evaluateDataFiles(evaluateFStatFiles);
}

void CGHandler::computeContactPoints() {
    //recompute contact point (necessary for old restart files)
    for (BaseInteraction* const c : getDPMBase()->interactionHandler)
    {
        if (c->getContactPoint().isNaN())
        {
            const Vec3D& p = c->getP()->getPosition();
            const Vec3D& i = c->getI()->getPosition();
            const BaseParticle* const PParticle = dynamic_cast<BaseParticle*>(c->getP());
            const BaseParticle* const IParticle = dynamic_cast<BaseParticle*>(c->getI());
            if (IParticle != nullptr)
            {
                const Vec3D branchVector = p - i;
                const Mdouble distance = branchVector.getLength();
                c->setNormal(branchVector / distance);
                c->setOverlap(PParticle->getRadius() + IParticle->getRadius() - distance);
                c->setDistance(distance);
                c->setContactPoint(p - (PParticle->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
            } else {
                const InfiniteWall* const IWall = dynamic_cast<InfiniteWall*>(c->getI());
                if (IWall != nullptr)
                {
                    const Mdouble dist = Vec3D::dot(i - p, IWall->getNormal());
                    const Mdouble r = dynamic_cast<BaseParticle*>(c->getP())->getRadius();
                    const Vec3D branch = (0.5 * (dist + r)) * IWall->getNormal();
                    c->setNormal(-IWall->getNormal());
                    c->setContactPoint(p + branch);
                } else
                {
                    const Vec3D IP = p - i;
                    const Mdouble distance = Vec3D::getLength(IP);
                    c->setNormal(IP / distance);
                    c->setContactPoint(i);
                }
            }
            //logger(INFO,"c%",c->getContactPoint());
            static bool firstTime = true;
            if (firstTime)
            {
                cgLogger(WARN,
                         "No contactPoint information and has to be recomputed; "
                         "for complex walls and non-spherical particles, midpoint rule is assumed; "
                         "periodic walls can create errors");
                firstTime = false;
            }
        }
    }
}

void CGHandler::evaluateRestartFiles()
{

#ifdef MERCURY_USE_MPI
    //Make sure that the number of processors is equal to the number of processors used for the run
    MPIContainer& communicator = MPIContainer::Instance();
    std::vector<unsigned> numberOfDomains = this->getDPMBase()->getNumberOfDomains();
    int numberOfRequiredProcessors = numberOfDomains[0]*numberOfDomains[1]*numberOfDomains[2];
    if(!(numberOfRequiredProcessors == communicator.getNumberOfProcessors()))
    {
        if (communicator.getProcessorID() == 0)
        {
            logger(ERROR,"Please re-run the program with % cores",numberOfRequiredProcessors);
        }
        else
        {
            std::exit(-1);
        }
    }
#endif
    // reset counters so reading begins with the first data/fstat file
    DPMBase* dpm = getDPMBase();
    //dpm->interactionHandler.clear();
    dpm->dataFile.setCounter(0);
    dpm->fStatFile.setCounter(0);

    initialise();

    // read data file
    do
    {
        //recompute contact point (necessary for old restart files)
        computeContactPoints();

        cgLogger(INFO, "Successfully read %, t=%, Np=%, Nc=%", dpm->restartFile.getFullName(), dpm->getTime(),
                 dpm->particleHandler.getSize(), dpm->interactionHandler.getNumberOfObjects());

        // Check if we still need to evaluate data files
        bool unfinished = true;
        for (BaseCG* it : *this)
        {
            //if we are below timeMax and the next time step should be written
            if (it->getTimeMax() < getDPMBase()->getTime())
            {
                std::cout << "cg time max: " << it->getTimeMax() << std::endl;
                std::cout << "dpm time max: " << dpm->getTime() << std::endl;
                unfinished = false;
            }
        }
        if (unfinished)
        {
            evaluate();
        } else
        {
            cgLogger(INFO, "Final time has been reached");
            break;
        }

        //the particle and wall handler is cleared here, because  BaseSpecies doesn't delete particles belonging to it
        dpm->particleHandler.clear();
        dpm->wallHandler.clear();

    } while (dpm->readRestartFile());

    finish();

    cgLogger(INFO, "Finished reading from %", dpm->dataFile.getFullName());
    dpm->dataFile.close();

}

void CGHandler::evaluateDataFiles(bool evaluateFStatFiles)
{

#ifdef MERCURY_USE_MPI
    //Make sure that the number of processors is equal to the number of processors used for the run
    MPIContainer& communicator = MPIContainer::Instance();
    std::vector<unsigned> numberOfDomains = this->getDPMBase()->getNumberOfDomains();
    int numberOfRequiredProcessors = numberOfDomains[0]*numberOfDomains[1]*numberOfDomains[2];
    if(!(numberOfRequiredProcessors == communicator.getNumberOfProcessors()))
    {
        if (communicator.getProcessorID() == 0)
        {
            logger(ERROR,"Please re-run the program with % cores",numberOfRequiredProcessors);
        }
        else
        {
            std::exit(-1);
        }
    }
#endif
    // reset counters so reading begins with the first data/fstat file
    DPMBase* dpm = getDPMBase();
    dpm->interactionHandler.clear();
    dpm->dataFile.setCounter(0);
    dpm->fStatFile.setCounter(0);

    initialise();

    // read data file
    while (dpm->readNextDataFile())
    {
        if (evaluateFStatFiles) dpm->readNextFStatFile();
        cgLogger(INFO, "Successfully read %, t=%, Np=%, Nc=%", dpm->dataFile.getFullName(), dpm->getTime(),
                 dpm->particleHandler.getSize(), dpm->interactionHandler.getNumberOfObjects());
        //cgLogger(INFO,"Successfully read %, t=%",dpm->fStatFile.getName(), dpm->getTime());

        // Check if we still need to evaluate data files
        bool unfinished = true;
        for (BaseCG* it : *this)
        {
            //if we are below timeMax and the next time step should be written
            if (it->getTimeMax() < getDPMBase()->getTime())
            {
                std::cout << "cg time max: " << it->getTimeMax() << std::endl;
                std::cout << "dpm time max: " << getDPMBase()->getTime() << std::endl;
                unfinished = false;
            }
        }
        if (unfinished)
        {
            evaluate();
        } else
        {
            cgLogger(INFO, "Final time has been reached");
            break;
        }

    }

    finish();

    cgLogger(INFO, "Finished reading from %", dpm->dataFile.getFullName());
    dpm->dataFile.close();

}
