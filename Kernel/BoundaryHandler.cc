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

#include "Math/Helpers.h"
#include "BoundaryHandler.h"
#include "Boundaries/BasePeriodicBoundary.h"
#include "Boundaries/AngledPeriodicBoundary.h"
#include "Boundaries/ChuteInsertionBoundary.h"
#include "Boundaries/CircularPeriodicBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/CubeDeletionBoundary.h"
#include "Boundaries/HopperInsertionBoundary.h"
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/ConstantMassFlowMaserBoundary.h"
#include "Boundaries/SubcriticalMaserBoundary.h"
#include "Boundaries/SubcriticalMaserBoundaryTEST.h"
#include "Boundaries/FluxBoundary.h"
#include "Boundaries/HeaterBoundary.h"
#include "Boundaries/BaseClusterInsertionBoundary.h"
#include "Boundaries/RandomClusterInsertionBoundary.h"
#include "Boundaries/FixedClusterInsertionBoundary.h"
#include "Boundaries/StressStrainControlBoundary.h"
#include "Boundaries/DropletBoundary.h"

///Constructor of the BoundaryHandler class. It creates and empty BoundaryHandler.
BoundaryHandler::BoundaryHandler()
{
    writeVTK_ = false;
    logger(DEBUG, "BoundaryHandler::BoundaryHandler() finished");
}

/*!
 * \param[in] BH The BoundaryHandler that has to be copied.
 * \details This is not a copy constructor! It just copies all BaseBoundary from
 *          the other handler into this handler, and clears all other variables.
 */
BoundaryHandler::BoundaryHandler(const BoundaryHandler& BH)
        : BaseHandler<BaseBoundary>()
{
    writeVTK_ = BH.writeVTK_;
    copyContentsFromOtherHandler(BH);
    logger(DEBUG, "BoundaryHandler::BoundaryHandler(const BoundaryHandler &BH) finished");
}

/*!
 * \param[in] rhs The BoundaryHandler on the right hand side of the assignment.
 * \details This is not a copy assignment operator! It just copies all BaseBoundary
 *          from the other handler into this handler, and clears all other variables.
 */
BoundaryHandler& BoundaryHandler::operator=(const BoundaryHandler& rhs)
{
    if (this != &rhs)
    {
        clear();
        copyContentsFromOtherHandler(rhs);
    }
    logger(DEBUG, "BoundaryHandler BoundaryHandler::operator =(const BoundaryHandler& rhs)");
    
    return *this;
}

///Default destructor. Note that the delete for all boundaries is done in the BaseHandler.
BoundaryHandler::~BoundaryHandler()
{
    logger(DEBUG, "BoundaryHandler::~BoundaryHandler() finished");
}

///\param[in] P A pointer to the BaseBoundary (or derived class) that has to be added.
///Add the object and tell the object that this is his handler.
void BoundaryHandler::addObject(BaseBoundary* P)
{
    //Puts the boundary in the Boundary list
    BaseHandler<BaseBoundary>::addObject(P);
    //set the handler pointer
    P->setHandler(this);

#ifdef MERCURY_USE_MPI
    //Different manner of treating periodic boundaries in MPI: hence there is a periodicBoundaryHandler
    BasePeriodicBoundary* upcast = dynamic_cast<BasePeriodicBoundary*>(P);
    if (upcast != nullptr)
    {
        getDPMBase()->periodicBoundaryHandler.addObject(upcast);
    }
#endif

}

BaseBoundary* BoundaryHandler::createObject(const std::string& type)
{
    //Note that compare returns 0 if the strings are the same.
    if (type == "AngledPeriodicBoundary")
    {
        return new AngledPeriodicBoundary;
    }
    else if (type == "ChuteInsertionBoundary")
    {
        return new ChuteInsertionBoundary;
    }
    else if (type == "CubeInsertionBoundary")
    {
        return new CubeInsertionBoundary;
    }
    else if (type == "CubeDeletionBoundary")
    {
        return new CubeDeletionBoundary;
    }
    else if (type == "CircularPeriodicBoundary")
    {
        return new CircularPeriodicBoundary;
    }
    else if (type == "DeletionBoundary")
    {
        return new DeletionBoundary;
    }
    else if (type == "HopperInsertionBoundary")
    {
        return new HopperInsertionBoundary;
    }
    else if (type == "PeriodicBoundary")
    {
        return new PeriodicBoundary;
    }
    else if (type == "ConstantMassFlowMaserBoundary" || type == "MaserBoundary")
    {
        return new ConstantMassFlowMaserBoundary;
    }
    else if (type == "SubcriticalMaserBoundary")
    {
        return new SubcriticalMaserBoundary;
    }
    else if (type == "SubcriticalMaserBoundaryTEST")
    {
        return new SubcriticalMaserBoundaryTEST;
    }
    else if (type == "LeesEdwardsBoundary")
    {
        return new LeesEdwardsBoundary;
    }
    else if (type == "MPIDomainBoundary")
    {
        return new PeriodicBoundary;
    }
    else if (type == "FluxBoundary")
    {
        return new FluxBoundary;
    }
    else if (type == "HeaterBoundary")
    {
        return new HeaterBoundary;
    }
    else if (type == "StressStrainControlBoundary")
    {
        return new StressStrainControlBoundary;
    }
    else if (type == "BaseClusterInsertionBoundary")
    {
        return new BaseClusterInsertionBoundary;
    }
    else if (type == "RandomClusterInsertionBoundary")
    {
        return new RandomClusterInsertionBoundary;
    }
    else if (type == "FixedClusterInsertionBoundary")
    {
        return new FixedClusterInsertionBoundary;
    }
    else if (type == "DropletBoundary")
    {
        return new DropletBoundary;
    }
    else if (type == "normal") //for backward compatibility (before svnversion ~2360)
    {
        return new PeriodicBoundary;
    }
    else
    {
        logger(WARN, "Boundary type: % not understood in restart file.", type);
        return nullptr;
    }
}

///\param[in] is The input stream from which the information is read.
///First read the type of boundary, then compare the type to all existing types. 
///When the correct type is found, read it with the >> operator, copy it and add it to the handler.
void BoundaryHandler::readAndAddObject(std::istream& is)
{
    //Note that compare returns 0 if the strings are the same.
    std::string type;
    is >> type;
    logger(VERBOSE, "BoundaryHandler::readAndAddObject(is): restarting a boundary of type %.", type);
    if (type == "normal")
    {
        readOldObject(is);
    }
    else
    {
        BaseBoundary* boundary = createObject(type);
        if (boundary == nullptr)
        {
            std::string line;
            getline(is, line);
            logger(WARN, "This boundary could not be read and is ignored:\n%%", type, line);
            return;
        }
        boundary->setHandler(this);
        is >> *boundary;
        addObject(boundary);
    }
}

///\param[in] is The input stream from which the information is read.
///Get the normal, position left and position right for a periodic boundary from the 
///stream is, and construct a new periodic boundary from it.
///The boundaries that are written like that are outdated, this function is there for backward compatability.
void BoundaryHandler::readOldObject(std::istream& is)
{
    //read in next line
    std::stringstream line;
    helpers::getLineFromStringStream(is, line);
    
    std::string dummy;
    Vec3D normal;
    Mdouble positionLeft, positionRight;
    
    PeriodicBoundary periodicBoundary;
    line >> normal >> dummy >> positionLeft >> dummy >> positionRight;
    periodicBoundary.set(normal, positionLeft, positionRight);
    copyAndAddObject(periodicBoundary);
}

/// \return The string "BoundaryHandler"
std::string BoundaryHandler::getName() const
{
    return "BoundaryHandler";
}

void BoundaryHandler::boundaryActionsBeforeTimeLoop()
{
    for (BaseBoundary* b : objects_)
    {
        b->actionsBeforeTimeLoop();
    }
}

