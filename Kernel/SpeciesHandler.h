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


#ifndef SPECIESHANDLER_H
#define SPECIESHANDLER_H

#include "BaseHandler.h"
#include "Species/ParticleSpecies.h"

/// \brief Container to store all ParticleSpecies
/// \details The SpeciesHandler is a container to store all ParticleSpecies. 
/// It is implemented by a vector of pointers to ParticleSpecies.
class SpeciesHandler final : public BaseHandler<ParticleSpecies>
{
public:
    /// \brief Default constructor, it creates an empty SpeciesHandler.
    SpeciesHandler();
    
    /// \brief Constructor that copies all species and the pointer to the DPMBase from the given SpeciesHandler.
    SpeciesHandler(const SpeciesHandler& other);
    
    /// \brief Assignment operator that copies all species and the pointer to the DPMBase from the given SpeciesHandler.
    SpeciesHandler& operator=(const SpeciesHandler& rhs);
    
    /// \brief Destructor, it destructs the SpeciesHandler and all ParticleSpecies it contains.
    ~SpeciesHandler() override;
    
    /// \brief Adds a new ParticleSpecies to the SpeciesHandler.
    void addObject(ParticleSpecies* S) override;
    
    void clear() override
    {
        BaseHandler<ParticleSpecies>::clear();
        mixedObjects_.clear();
    }
    
    ///\brief Remove the ParticleSpecies with given id.
    void removeObject(unsigned int index) override;
    
    /// \brief Reads Species data into the SpeciesHandler from restart file.
    void readAndAddObject(std::istream& is) override;
    
    /// \brief Reads ParticleSpecies into the SpeciesHandler from old-style restart data.
    ParticleSpecies* readOldObject(std::istream& is);
    
    /// \brief Gets the Id of the behaviour between two given species.
    unsigned int getMixedId(unsigned int id1, unsigned int id2) const;
    
    template<typename U>
    typename std::enable_if<!std::is_pointer<typename U::MixedSpeciesType>::value, typename U::MixedSpeciesType*>::type
    getMixedObject(const U* S, const U* T)
    {
        return static_cast<typename U::MixedSpeciesType*>(getMixedObject(S->getIndex(), T->getIndex()));
    }
    
    /// \brief Gets the mixed object that is constructed from two given species.
    BaseSpecies* getMixedObject(unsigned int id1, unsigned int id2);
    
    /// \brief Returns a pointer to the vector of all mixed objects.
    const std::vector<BaseSpecies*>& getMixedObjects() const;
    
    /// \brief Write all the species and mixed species to an output stream.
    virtual void write(std::ostream& os) const;
    
    /// \brief Returns the name of the handler, namely the string "SpeciesHandler".
    std::string getName() const override;
    
    /// \brief Check if angular DOF have to be used
    bool useAngularDOFs();

private:
    ///The list of pointers to the mixed species
    std::vector<BaseSpecies*> mixedObjects_;
};

#endif

