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

/*!
 * \file DomainHandler.h
 *
 * Class DomainHandler
 */
#ifndef DOMAINHANDLER_H
#define DOMAINHANDLER_H

#include "BaseHandler.h"
#include "Boundaries/BaseBoundary.h"
#include "Domain.h"

/*!
 * \class DomainHandler
 * \brief Container to store all Domain.
 * \details The DomainHandler is a container to store all Domain. 
 * It is implemented by a vector of pointers to domain. Additionally
 * it also contains an index that points to the current domain of the
 * processorID. Such that the processor knows which Domain it has to compute
 */
class DomainHandler final : public BaseHandler<Domain>
{
public:
    /*!
     * \brief Default constructor, it creates an empty DomainHandler.
     */
    DomainHandler();
    
    /*!
     * \brief Constructor that copies all Domain it contains.
     */
    DomainHandler(const DomainHandler& DH);
    
    /*!
     * \brief Assignment operator. 
     */
    DomainHandler& operator=(const DomainHandler& rhs);
    
    /*!
     * \brief Destructor, it destructs the DomainHandler and all Domain it contains.
     */
    ~DomainHandler() final;
    
    /*!
     * \brief Creates a Cartesian square mesh in 3D
     */
    void createMesh(std::vector<Mdouble>& simulationMin, std::vector<Mdouble>& simulationMax,
                    std::vector<unsigned>& numberOfDomains, bool open);
    
    /*!
     * \brief Recursive function that creates the domains for a 3D mesh
     */
    void createDomains(std::vector<Mdouble> domainMin, std::vector<Mdouble> domainMax,
                       std::vector<unsigned>& globalMeshIndex, std::vector<unsigned>& numberOfDomains, int dimCounter,
                       std::vector<Mdouble>& meshSize, bool open);
    
    /*!
     * \brief Adds a Domain to the DomainHandler. 
     */
    void addObject(Domain* D) final;
    
    /*!
     * \brief \todo Still has to be implemented 
     */
    void readAndAddObject(std::istream& is) final;
    
    /*!
     * \brief \todo Still has to be implemented 
     */
    void readOldObject(std::istream& is);
    
    /*!
     * \brief \todo Still has to be implemented 
     */
    std::string getName() const final;
    
    /*!
     * \brief This sets a domain to the processor
     */
    void setCurrentDomainIndex(unsigned int index);
    
    /*!
     * \brief Gets the domain assigned to the processor
     */
    Domain* getCurrentDomain();
    const Domain* getCurrentDomain() const;

    /*!
     * \brief Gets the domain index assigned to the processor
     */
    unsigned int getCurrentDomainIndex() const;
    
    /*!
     * \brief Sets the number of domains in the domain handler
     */
    void setNumberOfDomains(std::vector<unsigned>& numberOfdomains);
    
    /*!
     * \brief  Gets the number of domains in the domain handler
     */
    std::vector<unsigned> getNumberOfDomains();
    
    /*!
     * \brief Sets the interaction distance of the domain handler
     */
    void setInteractionDistance(Mdouble interactionDistance);
    
    /*!
     * \brief Gets the interaction distance of the domain handler
     */
    Mdouble getInteractionDistance();
    
    /// \todo MX: function under construction
    ///\todo TW@Marnix should this be unsigned int?
    int getParticleDomainGlobalIndex(BaseParticle* particle);
    
    ///\todo TW@Marnix should this be unsigned int?
    int getParticleProcessor(int globalIndex);
    
    ///\todo TW@Marnix should this be unsigned int?
    Domain* getParticleDomain(int globalIndex);
    
    void updateStatus(std::set<BaseParticle*>& particlesToBeDeleted);
    
    void updateVelocity();
    
    void addNewParticles();
    
    void initialise();

private:
    
    /*!
     * \brief Index to the particular domain of this process
     */
    unsigned int currentDomainIndex_;
    
    /*!
     * \brief look-up table to find the processor of a domain given the globalIndex of the domain
     */
    std::vector<int> globalIndexToProcessorList_;
    
    /*!
     * \brief vector containing the number of domains in Cartesian direction
     */
    std::vector<unsigned> numberOfDomains_;
    
    /*!
     * \brief Interaction distance between a domain boundary and the communication zone boundary
     */
    Mdouble interactionDistance_;
};

#endif

