//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include "MpiContainer.h"
#include "Logger.h"

#ifdef MERCURYDPM_USE_MPI
#include <mpi.h>
#include "Interactions/BaseInteraction.h"
#include "Interactions/Interaction.h"
#include "Particles/BaseParticle.h"
#include "MpiDataClass.h"
#endif

#include "SpeciesHandler.h"
#include "GeneralDefine.h"

/*! 
 * \brief Initialise the communicator after MPI has been initalised.
 */
MPIContainer::MPIContainer()
{
#ifdef MERCURYDPM_USE_MPI
    int Mpi_init_flag = 0;
    MPI_Initialized(&Mpi_init_flag);
    if(!Mpi_init_flag)
    {
        logger(FATAL,"MPI should be initialised before calling the MPIContainer constructor");
    }
    //A personal communicator will be created to ensure we don't meddle with communicators of other libraries
    MPI_Group groupID;
    MPI_Comm_group(MPI_COMM_WORLD, &groupID);
    MPI_Comm_create(MPI_COMM_WORLD,groupID,&communicator_);
    MPI_Comm_rank(MPI_COMM_WORLD,&processorID_);
    MPI_Comm_size(MPI_COMM_WORLD,&numberOfProcessors_);

#else
    numberOfProcessors_ = 1;
    processorID_ = 0;
#endif
}

/*
 * \brief Initialises MercuryMPIType MPITypes 
 * \details This function intialises five different Mercury data types that have to
 * be send over the MPI interface. Most of the data types are just structs containing the raw
 * required data. The interaction data type however is a very complex datatype
 * as it can vary over every simulation what type of interaction is set.
 * Currently it is only possible to have one type of interaction in the MPI (and that should be sufficient)
 * \param[in] speciesHandler Handles all the species, required to create a prototype interaction
 */
void MPIContainer::initialiseMercuryMPITypes(const SpeciesHandler& speciesHandler)
{
#ifdef MERCURYDPM_USE_MPI
    //Note: Important that the MPI type creation is done in the order given by the enum
    MPIParticle dummyParticle;
    MPIParticlePosition dummyPosition;
    MPIParticleVelocity dummyVelocity;
    MPIParticleForce dummyForce;
    createMercuryMPIType(dummyParticle, MercuryMPIType::PARTICLE);
    createMercuryMPIType(dummyPosition, MercuryMPIType::POSITION);
    createMercuryMPIType(dummyVelocity, MercuryMPIType::VELOCITY);
    createMercuryMPIType(dummyForce,    MercuryMPIType::FORCE);
    
    //Obtain the correct history force class
    //NOTE: only works for one type of interaction
    if (dataTypes_.size() == 4)
    {
        //Create a dummy interaction to get a grip on the size of the MPI interaction class
        const BaseSpecies* species = speciesHandler.getObject(0);
        BaseInteraction* emptyInteraction = species->getEmptyInteraction();
        emptyInteraction->createMPIType();
        species->deleteEmptyInteraction(emptyInteraction);
    }
#endif
}

/*!
 * \brief Obtain the number of processors
 * \return Returns the number of processors
 */
std::size_t MPIContainer::getNumberOfProcessors() const
{
    return numberOfProcessors_;
}

/*!
 * \brief Obtain the current processor ID
 * \return Returns the processor ID in the communication group
 */
std::size_t MPIContainer::getProcessorID()
{
    return processorID_;
}


#ifdef MERCURYDPM_USE_MPI

/*!
 * \brief Get the communicator used for MPI commands.
 * \details The MPI::Intracomm is required for parsing information
 * between processors
 * \return Returns the private MercuryDPM MPI communicator
 */
MPI_Comm& MPIContainer::getComm()
{
    return communicator_;
}

#endif

/*!
 * \brief Inialises the MPI library
 */
void initialiseMPI()
{
#ifdef MERCURYDPM_USE_MPI
    //Check if MPI is already initialised
    int Mpi_init_flag = 0;
    MPI_Initialized(&Mpi_init_flag);
    if(!Mpi_init_flag)
    {
        MPI_Init(NULL,NULL);
        MPIContainer& communicator = MPIContainer::Instance();
        if (PROCESSOR_ID  == 0)
        {
            std::cout << "MPI has been initialised" << std::endl;
        }

        //MPI should be finalised at the end of any program.
        std::atexit([]()
        {
            MPIContainer& communicator = MPIContainer::Instance();
            communicator.deleteMercuryMPITypes();
            MPI_Finalize();
            if (PROCESSOR_ID == 0)
            {
                std::cout << "MPI has been finalised" << std::endl;
            }
        });
    }
#endif
}

