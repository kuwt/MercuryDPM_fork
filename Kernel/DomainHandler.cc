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


#include <iostream>
#include <string>
#include <cstdlib>
#include <Math/Helpers.h>
#include "DomainHandler.h"
#include "Boundaries/BaseBoundary.h"
#include "DPMBase.h"
#include <limits>

/*!
 * \details Constructor of the DomainHandler class. It creates and empty 
 *          DomainHandler.
 */
DomainHandler::DomainHandler()
{
    currentDomainIndex_ = 0;
    interactionDistance_ = 0.0;
    logger(DEBUG, "DomainHandler::DomainHandler() finished");
}

/*!
 * \param[in] DH The DomainHandler that has to be copied. 
 * \details This is not a copy constructor! It copies the DPMBase and all 
 *          Domain, and sets the other variables to 0.
 */
DomainHandler::DomainHandler(const DomainHandler& DH)
        : BaseHandler<Domain>()
{
    copyContentsFromOtherHandler(DH);
    logger(DEBUG, "DomainHandler::DomainHandler(const DomainHandler &DH) finished");
}

/*!
 * \param[in] rhs The DomainHandler on the right hand side of the assignment.
 * \details This is not a copy assignment operator! It copies the DPMBase and all 
 *          Domain, and sets the other variables to 0.
 */
DomainHandler& DomainHandler::operator=(const DomainHandler& rhs)
{
    if (this != &rhs)
    {
        clear();
        copyContentsFromOtherHandler(rhs);
    }
    logger(DEBUG, "DomainHandler DomainHandler::operator =(const DomainHandler& rhs)");
    return *this;
}

/*!
 * \details All Domain are destroyed, the DomainHandler afterwards.
 */
DomainHandler::~DomainHandler()
{
    logger(DEBUG, "DomainHandler::~DomainHandler() finished");
}

/*!
 * \brief Creates a cartesian square mesh in 3D
 * \details By default the number of domains should be 1 in each direction. The mesh creates the required 3D domains
 * and then disables any domain boundaries that are not connected to another domain. Lookup tables are created between domains
 * and ranks.
 * \param[in] simulationMin The minimum values of the simulation domain in cartesian coordinates
 * \param[in] simulationMax The maximum values of the simulation domain in cartesian coordinates
 * \param[in] numberOfDomains The number of domains in cartesian direction
 * \param[in] open Determines if the boundary domains have infinite limits (true) or not (false)
 */
void DomainHandler::createMesh(std::vector<Mdouble>& simulationMin, std::vector<Mdouble>& simulationMax,
                               std::vector<unsigned>& numberOfDomains, bool open)
{
    //Clear the objects in the list, we're gonna make a new mesh, owyeah
    objects_.clear();
    
    //Create the mesh
    setNumberOfDomains(numberOfDomains);
    
    std::vector<unsigned> globalMeshIndex(3);
    //Recursive function creating the domains
    int dimCounter = 3;
    //Compute the mesh size in each direction
    std::vector<Mdouble> meshSize;
    for (int d = 0; d < 3; d++)
    {
        meshSize.push_back((simulationMax[d] - simulationMin[d]) / numberOfDomains[d]);
    }
    //Create the domains
    createDomains(simulationMin, simulationMax, globalMeshIndex, numberOfDomains, dimCounter, meshSize, open);
    
    //Create lookup tables
    for (Domain* domain : objects_)
    {
        domain->setRank(domain->getGlobalIndex());
        domain->createLookUpTable();
    }
    
    //Disable boundaries that dont require communication
    for (Domain* domain : objects_)
    {
        domain->disableBoundaries();
    }
    
    //Output the result
    std::string meshType;
    if (open)
    {
        meshType = "open";
    }
    else
    {
        meshType = "closed";
    }
    if (PROCESSOR_ID == 0)
    {
        logger(INFO, "A simulation mesh has been created with % number of domains and % boundaries.", this->getSize(),
               meshType);
    }
}

/*!
 * \brief Recursive function that creates the domains for a 3D mesh.
 * \details For a given dimension, dimCounter, compute the size of the domain, and accordingly
 * create domains with the correct minimum and maximum value in the dimCounter-direction.
 * When all bounds of a domain are set, i.e. dimCounter is zero, the domain is actually created.
 * \param[in,out] domainMin The minimum simulation values in Cartesian direction as input and the minimum domain values as output/input
 * \param[in] domainMax The maximum simulation values in Cartesian direction as input and the maximum domain values as output/intput
 * \param[in] globalMeshIndex A Vector containing the mesh indices, i.e. (i,j,k) = (3,5,2)
 * \param[in] numberOfDomains A vector containg the number of domoains in Cartesian direction
 * \param[in] dimCounter The dimension counter that keeps track in whicn direction the rescursive function is working
 * \param[in] open determines if a domain boundary is open (inf) or closed
 */
void DomainHandler::createDomains(std::vector<Mdouble> domainMin, std::vector<Mdouble> domainMax,
                                  std::vector<unsigned>& globalMeshIndex, std::vector<unsigned>& numberOfDomains,
                                  int dimCounter, std::vector<Mdouble>& meshSize, bool open)
{
    //
    if (dimCounter == 0)
    {
        //Create a new domain
        Domain domain(globalMeshIndex);
        domain.setHandler(this);
        const std::vector<double>& domainBoundMin = domainMin;
        const std::vector<double>& domainBoundMax = domainMax;
        domain.setBounds(domainBoundMin, domainBoundMax, true);
        this->copyAndAddObject(domain);
    }
    else
    {
        dimCounter--;
        //Compute size of a domain in the dimCounter'th direction.
        Mdouble boundLeft = domainMin[dimCounter];
        //Starting with the bound left, create the number of domains in the given dimCounter-direction
        for (int i = 0; i < numberOfDomains[dimCounter]; i++)
        {
            globalMeshIndex[dimCounter] = i;
            domainMin[dimCounter] = boundLeft + i * meshSize[dimCounter];
            domainMax[dimCounter] = boundLeft + (i + 1) * meshSize[dimCounter];
            if ((i == 0) && (open))
            {
                domainMin[dimCounter] = -constants::inf;
            }
            
            if ((i == numberOfDomains[dimCounter] - 1) && (open))
            {
                domainMax[dimCounter] = constants::inf;
            }
            
            //Start recursively a new createDomains function
            createDomains(domainMin, domainMax, globalMeshIndex, numberOfDomains, dimCounter, meshSize, open);
        }
    }
}

/*!
 * \param[in] D A pointer to the Domain that has to be added. 
 * \details Adds the object to the DomainHandler and sets the DomainHandler
 * pointer in the Domain to this DomainHandler.
 */
void DomainHandler::addObject(Domain* D)
{
    //Puts the particle in the Particle list
    BaseHandler<Domain>::addObject(D);
    //set the particleHandler pointer
    D->setHandler(this);
}

/*!
 * \brief reads a domain object 
 * \details There is no need to read domain object, they can be computed easily
 */
void DomainHandler::readAndAddObject(std::istream& is)
{
}

/*!
 * \brief reads an old domain object 
 * \details There is no need to read a domain object, they can be computed easily
 */
void DomainHandler::readOldObject(std::istream& is)
{
}

/*!
 * \brief returns the name of the class 
 */
std::string DomainHandler::getName() const
{
    return "DomainHandler";
}

/*!
 * \param[in] index An integer to a Domain in the DomainHandler 
 * \details Sets an index to a Domain in the DomainHandler which belongs to 
 * the current processor.
 */
void DomainHandler::setCurrentDomainIndex(unsigned int index)
{
    currentDomainIndex_ = index;
}

/*!
 * \brief Returns a pointer to the active domain
 * \details The active domain is the domain that was assigned to the processor 
 * \return Pointer to a Domain on which the processor has to do computations
 */
Domain* DomainHandler::getCurrentDomain()
{
    return getObject(currentDomainIndex_);
}

const Domain* DomainHandler::getCurrentDomain() const
{
    return getObject(currentDomainIndex_);
}

/*!
 * \details Gets the Domain Index in the vector of the DomainHandler 
 * assigned to this processor 
 * \return Index of a Domain on which the processor has to do computations
 */
unsigned int DomainHandler::getCurrentDomainIndex() const
{
    return currentDomainIndex_;
}

/*!
 * \brief Sets the number of domains in the domain handler
 * \details The number of domains is given as a vector in the form of (nx,ny,nz)
 * where nx,ny,nz are the number of domains in their respective direction
 * \param[in} numberOfDomains the number of domains in Cartesian direction
 */
void DomainHandler::setNumberOfDomains(std::vector<unsigned>& numberOfDomains)
{
    numberOfDomains_ = numberOfDomains;
}

/*!
 * \brief Gets the number of domains in the domain handler
 * \details The number of domains is given as a vector in the form of (nx,ny,nz)
 * where nx,ny,nz are the number of domains in their respective direction
 * \return The number of domains in Cartesian direction
 */
std::vector<unsigned> DomainHandler::getNumberOfDomains()
{
    return numberOfDomains_;
}

/*!
 * \brief Sets the interaction distance of the domainHandler
 * \details The interaction distance is used to define how large the
 * communication zone of the domains have to be. generally twice the size of the
 * largest interaction radius
 * \param[in] interactionDistance Interaction distance between domain boundary and communication zone boundary
 */
void DomainHandler::setInteractionDistance(Mdouble interactionDistance)
{
    //Update the interaction distance
    interactionDistance_ = interactionDistance;
    
    //Check if the domainsize is not too small
    Domain* domain = getCurrentDomain();
    logger.assert_always((domain->getDomainMax()[0] - domain->getDomainMin()[0]) > 2 * interactionDistance_,
                         "Size of the domain in x-direction is smaller than communication zone. Size: %, communication zone: %",
                         (domain->getDomainMax()[0] - domain->getDomainMin()[0]), 2 * interactionDistance_);
    logger.assert_always((domain->getDomainMax()[1] - domain->getDomainMin()[1]) > 2 * interactionDistance_,
                         "Size of the domain in y-direction is smaller than communication zone. Size: %, communication zone: %",
                         (domain->getDomainMax()[1] - domain->getDomainMin()[1]), 2 * interactionDistance);
    logger.assert_always((domain->getDomainMax()[2] - domain->getDomainMin()[2]) > 2 * interactionDistance_,
                         "Size of the domain in z-direction is smaller than communication zone. Size: %, communication zone: %",
                         (domain->getDomainMax()[2] - domain->getDomainMin()[2]), 2 * interactionDistance_);
}

/*!
 * \brief Gets the interaction distance of the communication zone
 * \details The interaction distance is used to define how large the
 * communication zone of the domains have to be. generally twice the size of the
 * largest interaction radius
 * \return Returns the interaction distance of the communication zone
 */
Mdouble DomainHandler::getInteractionDistance()
{
    return interactionDistance_;
}

/// \todo MX: This function is still under development the goal of this function is to obtain the globalIndex of the domain the particle is located in
int DomainHandler::getParticleDomainGlobalIndex(BaseParticle* particle)
{
    //Step 1: obtain values i,j,k by looking at the position
    //TODO this could possibly be stored in the domainHandler to save computational power
    std::vector<int> decompositionVector(3);
    
    int i, j, k;
    Mdouble dx = (getDPMBase()->getXMax() - getDPMBase()->getXMin()) / getDPMBase()->getNumberOfDomains()[0];
    Mdouble dy = (getDPMBase()->getYMax() - getDPMBase()->getYMin()) / getDPMBase()->getNumberOfDomains()[1];
    Mdouble dz = (getDPMBase()->getZMax() - getDPMBase()->getZMin()) / getDPMBase()->getNumberOfDomains()[2];
    
    //x-direction
    if ((particle->getPosition().X - getDPMBase()->getXMin()) <= 0)
    {
        i = 0;
    }
    else
    {
        //Compute the relative domain cell it is in
        i = floor((particle->getPosition().X - getDPMBase()->getXMin()) / dx);
        
        //Make sure we dont go over our number of domain bounds
        if (i >= getDPMBase()->getNumberOfDomains()[0] - 1)
        {
            i = getDPMBase()->getNumberOfDomains()[0] - 1;
        }
        if (i < 0)
        {
            i = 0;
        }
    }
    
    //y-direction
    if ((particle->getPosition().Y - getDPMBase()->getYMin()) <= 0)
    {
        j = 0;
    }
    else
    {
        j = floor((particle->getPosition().Y - getDPMBase()->getYMin()) / dy);
        
        //Make sure we dont go over our number of domain bounds
        if (j >= getDPMBase()->getNumberOfDomains()[1] - 1)
        {
            j = getDPMBase()->getNumberOfDomains()[1] - 1;
        }
        if (i < 0)
        {
            j = 0;
        }
    }
    
    //z-direction
    if ((particle->getPosition().Z - getDPMBase()->getZMin()) <= 0)
    {
        k = 0;
    }
    else
    {
        k = floor((particle->getPosition().Z - getDPMBase()->getZMin()) / dz);
        //Make sure we dont go over our number of domain bounds
        if (k >= getDPMBase()->getNumberOfDomains()[2] - 1)
        {
            k = getDPMBase()->getNumberOfDomains()[2] - 1;
        }
        if (k < 0)
        {
            k = 0;
        }
        
    }
    
    //Step 2: obtain the processor number
    int globalIndex = i +
                      getDPMBase()->getNumberOfDomains()[0] * j +
                      getDPMBase()->getNumberOfDomains()[1] * getDPMBase()->getNumberOfDomains()[0] * k;
    return globalIndex;
}

int DomainHandler::getParticleProcessor(int globalIndex)
{
    return globalIndex;//ToProcessorList_[globalIndex];
}

Domain* DomainHandler::getParticleDomain(int globalIndex)
{
    return getObject(globalIndex);
}

void DomainHandler::updateStatus(std::set<BaseParticle*>& particlesToBeDeleted)
{
    getCurrentDomain()->updateStatus(particlesToBeDeleted);
}

void DomainHandler::updateVelocity()
{
    getCurrentDomain()->updateVelocity();
}

void DomainHandler::addNewParticles()
{
    getCurrentDomain()->addNewParticles();
}

void DomainHandler::initialise()
{
    //Create a single domain
    Domain domain;
    domain.setHandler(this);
    std::vector<double> domainBoundMin = {-constants::inf, -constants::inf, -constants::inf};
    std::vector<double> domainBoundMax = {constants::inf, constants::inf, constants::inf};
    domain.setBounds(domainBoundMin, domainBoundMax, false);
    this->copyAndAddObject(domain);
    currentDomainIndex_ = 0;
}
