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

/*!
 * \file MpiContainer.h
 *
 * Class MPIContainer
 */
#ifndef MPICONTAINER_H_
#define MPICONTAINER_H_

#include <cstddef>
#include <functional>
#include <vector>
#include "Math/Vector.h"

#ifdef MERCURYDPM_USE_MPI
#include <mpi.h>
#endif

#ifdef MERCURYDPM_FORCE_ASSERTS
#define MERCURYDPM_ASSERTS true
#else
#ifdef MERCURYDPM_NO_ASSERTS
#define MERCURYDPM_ASSERTS false
#else
#ifdef NDEBUG
#define MERCURYDPM_ASSERTS false
#else
#define MERCURYDPM_ASSERTS true
#endif
#endif
#endif

class SpeciesHandler;

/*!
 * \brief An enum that indicates what type of data is being send over MPI
 * \details MPI requires knowledge on the memory lay out of a data object that is being send.
 * Various types of data are being send in the parallel code, i.e. a whole particle or only the position.
 * This enum indicates what type is being send 
 */
enum MercuryMPIType
{
    PARTICLE = 0, POSITION = 1, VELOCITY = 2, FORCE = 3, INTERACTION = 4, SUPERQUADRIC = 5
};

/*!
 * \brief An enum that facilitates the creation of unique communication tags in the parallel code
 * \details The MercuryMPITag is the last digit of a communication tag in the parallel code. This ensures
 * that when various data types are queued simultationously, there is no communication confusion. Additionally
 * this is a useful feature for developers so they can trace back what message was being send if something went wrong. 
 */
enum MercuryMPITag
{
    PARTICLE_COUNT = 0,
    PARTICLE_DATA = 1,
    POSITION_DATA = 2,
    PERIODIC_POSITION_DATA = 3,
    VELOCITY_DATA = 4,
    INTERACTION_COUNT = 5,
    INTERACTION_DATA = 6,
    PERIODIC_COMPLEXITY = 7,
    PARTICLE_INDEX = 8,
    SUPERQUADRIC_DATA = 9
};

namespace Detail
{
#ifdef MERCURYDPM_USE_MPI
//convert integral data to the corresponding MPI type
template<typename T>
typename std::enable_if<std::is_integral<T>::value, MPI_Datatype>::type
toMPIType(T t)
{
    MPI_Datatype type;
    MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(T), &type);
    return type;
}

//convert floating point data to the corresponding MPI type
template<typename T>
typename std::enable_if<std::is_floating_point<T>::value, MPI_Datatype>::type
toMPIType(T t)
{
    MPI_Datatype type;
    MPI_Type_match_size(MPI_TYPECLASS_REAL,sizeof(T),&type);
    return type;
}
#endif
} /*detail*/

/*!
 * \brief Inialises the MPI library
 */
void initialiseMPI();

/*!
 * \class MPIContainer
 * \brief This class contains all information and functions required for communication between processors
 * \details This class contains of a newly created MPI_COMM communicator used communicate between processors.
 * The reason that MPI_COMM_WORLD is not used, is that other libraries might also be communicating in parallel
 * and you don't want to interfere with these communications. 
 * It is a signleton pattern, restricting to one instance of the class.
 * The class furthermore keeps track of the number of processors, processorID's and of asynchronous  send/receive requests. 
 */
class MPIContainer final
{
public:
    /// \brief fetch the instance to be used for communication
    /// \return The only instance of this class
    static MPIContainer& Instance()
    {
        static MPIContainer theInstance;
        return theInstance;
    }
    
    /*!
     * \brief Creates the MPI types required for communication of Mercury data through the MPI interface
     */
    void initialiseMercuryMPITypes(const SpeciesHandler& speciesHandler);
    
    /*! 
     * \brief Process all pending asynchronous communication requests before continuing
     * \details Each processor can commit a send and receive request at any time when is is 
     * finished with computing. At some point these requests have to be resolved in order to
     * keep order in the simulation. After all requests are resolved, a barrier ensures that
     * all processors wait untill all requests are resolved.
     */
    void sync()
    {
#ifdef MERCURYDPM_USE_MPI
        MPI_Waitall(pending_.size(),pending_.data(),MPI_STATUSES_IGNORE);
        pending_.clear();
        MPI_Barrier(communicator_);
#endif
    }
    
    /*! 
     * \brief Asynchronously send a scalar to some other processor.
     * \param[in] t the data
     * \param[in] to the processor to recieve the information
     * \param[in] tag an identifier for this specific send request. This must be unique among all send requests between
     * the previous synchronisation step and the next one. Exactly one recieve request must also provide this tag and
     * it must be done on the process specified by the 'to' parameter
     */
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    send(T& t, int to, int tag)
    {
#if MERCURYDPM_ASSERTS
        if (to == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Sending data to self!");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Request request;
        MPI_Isend(&t, 1, Detail::toMPIType(t), to, tag, communicator_, &request);
        pending_.push_back(request);

#endif
    }

/// \todo MX: type documentation. This one is used to send vectors of scalars across (hence the *t) 
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    send(T* t, int count, int to, int tag)
    {
#if MERCURYDPM_ASSERTS
        if (to == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Sending data to self!");
        }
        
        if (count == 0)
        {
            logger(WARN, "[MPI ERROR]: Sending zero data");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Request request;
        MPI_Isend(t, count, Detail::toMPIType(*t), to, tag, communicator_, &request);
        pending_.push_back(request);

#endif
    }
    
    /*! 
     * \brief asynchronously receive a scalar from some other processor.
     * \param[in,out] t the data
     * \param[in] from the processor that sends the information
     * \param[in] tag an identifier for this specific receive request. This must be unique among all receive requests between
     * the previous synchronisation step and the next one. Exactly one send request must also provide this tag and
     * it must be done on the process specified by the 'from' parameter
     */
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    receive(T& t, int from, int tag)
    {
#if MERCURYDPM_ASSERTS
        if (from == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Receiving data from self!");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Request request;
        MPI_Irecv(&t, 1, Detail::toMPIType(t), from, tag, communicator_, &request);
        pending_.push_back(request);

#endif
    }

/// \todo MX: type documentation. this is used to receive vectors of scalars accross   
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    receive(T* t, int count, int from, int tag)
    {
#if MERCURYDPM_ASSERTS
        if (from == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Receiving data fromself!");
        }
        
        if (count == 0)
        {
            logger(WARN, "[MPI ERROR]: Receiving zero data");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Request request;
        MPI_Irecv(&t, count, Detail::toMPIType(*t), from, tag, communicator_, &request);
        pending_.push_back(request);

#endif
    }
    
    
    /*!
     * \brief asynchronously send a list of MercuryMPITypes objects to some other processor.
     * \param[in,out] t the data, list of MPIType objects
     * \param[in] type the MPIType that is being send, for instance an MPIParticle
     * \param[in] count the number of objects that are being send
     * \param[in] to the processor that the data is send to
     * \param[in] tag a unique identifier that corresponds with a receive command by the receiving processor 
     */
    template<typename T>
    void send(T* t, MercuryMPIType type, int count, int to, int tag)
    {
        //std::cout << "[Process: " << processorID_ << "]: QUEUING SEND REQUEST with tag: " << tag << std::endl;
#if MERCURYDPM_ASSERTS
        if (to == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Sending data to self!");
        }
        
        if (count == 0)
        {
            logger(WARN, "[MPI ERROR]: Sending zero data");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Request request;
        MPI_Isend(t, count, dataTypes_[type], to, tag, communicator_, &request);
        pending_.push_back(request);


#endif
    }
    
    /*!
     * \brief asynchronously receive a list of MercuryMPIType objects from some other processor.
     * \param[in,out] t the data, list of MercuryMPIType objects
     * \param[in] type the MPIType that is being received, for instance an MPIParticle
     * \param[in] count the number of objects that are being send
     * \param[in] from the processor that sends the information
     * \param[in] tag a unique identifier that corresponds with a send command by the sending processor
     */
    template<typename T>
    void receive(T* t, MercuryMPIType type, int count, int from, int tag)
    {
        //std::cout << "[Process: " << processorID_ << "]: QUEUING RECEIVE REQUEST with tag: " << tag << std::endl;
#if MERCURYDPM_ASSERTS
        if (from == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Receiving data to self!");
        }
        
        if (count == 0)
        {
            logger(WARN, "[MPI ERROR]: Receiving zero data");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Request request;
        MPI_Irecv(t, count, dataTypes_[type], from, tag, communicator_, &request);
        pending_.push_back(request);

#endif
    }
    
    /*!
     * \brief synchronously send a list of scalars to another processor. 
     * the data should be received directly or the program will stall
     * \param[in,out] t the data, list of scalars
     * \param[in] count the number of scalars to be send
     * \param[in] to the processor that the data is being send to
     * \param[in] tag a unique identifier that corresponds with a receive command by the receiving processor
     */
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    directSend(T& t, int count, int to, int tag)
    {
#if MERCURYDPM_ASSERTS
        if (to == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Sending data to self!");
        }
        
        if (count == 0)
        {
            logger(WARN, "[MPI ERROR]: Sending zero data");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Ssend(&t, count, Detail::toMPIType(t), to, tag, communicator_);
#endif
    }
    
    
    template<typename T>
    void directSend(T* t, MercuryMPIType type, int count, int to, int tag)
    {
#if MERCURYDPM_ASSERTS
        if (to == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Sending data to self!");
        }
        
        if (count == 0)
        {
            logger(WARN, "[MPI ERROR]: Sending zero data");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Ssend(t,count,dataTypes_[type], to, tag, communicator_);
#endif
    }
    
    /*!
     * \brief synchronously receive a list of scalars from another processor.
     * if the send command has not been issued, this function will stall the program
     * \param[in,out] t the data, list of scalars
     * \param[in] count the number of scalars to be send
     * \param[in] from the processor that sends the information
     * \param[in] tag a unique identifier that corresponds with a send command by the sending processor
     */
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    directReceive(T& t, int count, int from, int tag)
    {
#if MERCURYDPM_ASSERTS
        if (from == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Receiving data from self!");
        }
        
        if (count == 0)
        {
            logger(WARN, "[MPI ERROR]: Receiving zero data");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Recv(&t, count, Detail::toMPIType(t), from, tag,communicator_, MPI_STATUS_IGNORE);
#endif
    }
    
    template<typename T>
    void directReceive(T* t, MercuryMPIType type, int count, int from, int tag)
    {
#if MERCURYDPM_ASSERTS
        if (from == processorID_)
        {
            logger(FATAL, "[MPI FATAL]: Receiving data to self!");
        }
        
        if (count == 0)
        {
            logger(WARN, "[MPI ERROR]: Receiving zero data");
        }
#endif
#ifdef MERCURYDPM_USE_MPI
        MPI_Recv(t, count, dataTypes_[type], from, tag, communicator_, MPI_STATUS_IGNORE);
#endif
    }
    
    /*!
     * \brief Gathers a scaler from all processors to a vector of scalars on the root
     * \details When a single processor needs to know a certain value on all processors, this
     * function will gathers them together at the root processor in an array of the size of the communicator size
     * \param[in,out] send_t the data that is being send, scalar value
     * \param[in,out] receive_t the data that is being received by the root, an array of scalars
     */
    template<typename T>
    void gather(T& send_t, T* receive_t)
    {
#ifdef MERCURYDPM_USE_MPI
        MPI_Gather(&send_t, 1, Detail::toMPIType(send_t), receive_t, 1, Detail::toMPIType(send_t), 0, communicator_);
#endif
    }
    
    /*!
     * \brief Broadcasts a scalar from the root to all other processors
     * \param[in,out] t scalar data that is being send by the root and received by the processors
     */
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    broadcast(T& t, int fromProcessor = 0)
    {
#ifdef MERCURYDPM_USE_MPI
        MPI_Bcast(&t,1,Detail::toMPIType(t),fromProcessor,communicator_);
#endif
    }
    
    /*!
     * \brief Broadcasts a scalar from the root to all other processors
     * \param[in,out] t scalar data that is being send by the root and received by the processors
     */
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    broadcast(T* t, int size, int fromProcessor)
    {
#ifdef MERCURYDPM_USE_MPI
        MPI_Bcast((void *)t,size,Detail::toMPIType(t[0]),fromProcessor,communicator_);

#endif
    }
    
    /*!
     * \brief Broadcasts an MercuryMPIType to all other processors
     * \param[in,out] t MercuryMPIType data that is being send by the root and received by the processors
     */
    template<typename T>
    void broadcast(T* t, MercuryMPIType type, int fromProcessor = 0)
    {
#ifdef MERCURYDPM_USE_MPI
        MPI_Bcast((void *)t,1,dataTypes_[type],fromProcessor,communicator_);
#endif
    }
    
    /*!
     * \brief Reduces a scalar on all processors to one scalar on a target processor
     * \details A scalar defined on all processors is reduced to one number by an operation
     * examples of operations are MPI::MAX, MPI::MIN and MPI::SUM, giving the maximum, minumum
     * or the summation of the given scalars. The resulting reduced scalar is then sent to
     * the target processor id, generally this is the root processor 0.
     * \param[in,out] t Scalar that needs to be collected and reduced to one scalar
     * \param[in] operation Operation that is performed on the collected scalars
     * \param[in] id Optional input, receiving processor of the reduced scalar
     */
#ifdef MERCURYDPM_USE_MPI
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    reduce(T& t, MPI_Op operation, int id = 0)
    {

        if(id == getProcessorID())
        {
            MPI_Reduce(MPI_IN_PLACE, &t, 1, Detail::toMPIType(t), operation, id, communicator_);
        }
        else
        {
            MPI_Reduce(&t, nullptr, 1, Detail::toMPIType(t), operation, id, communicator_);
        }
    }
#endif
    
    /*!
     * \brief AllReduces a scalar on all processors by a given MPI operation
     * \details A local scalar on all processors is reduced to one scalar as output
     * the reduction follows the operation rule given.
     * \param[in] send_t the scalar that is send to be reduced
     * \param[out] receive_t the reduced scalar 
     * \param[in] operation The operation that is performed on all local numbers to reduce it to a single number
     */
#ifdef MERCURYDPM_USE_MPI
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    allReduce(T& send_t, T& receive_t, MPI_Op operation)
    {
        MPI_Allreduce(&send_t, &receive_t, 1, Detail::toMPIType(send_t), operation,communicator_);  
    }
#endif
    
    /*!
     * \brief allGather takes a (different) scalar from all processors and returns a vector with all scalars
     * \details sometimes it is required to share a local scalar such as number of particles to all other processors.
     * With allGather all processors now now the local scalar of all other processors
     * \param[in] send_t The scalar that is send to all other processors
     * \param[in] send_count the number of scalars send to other processors
     * \param[in,out] receive_t A vector of scalars that contain all other scalars from other processors
     * \param[in] receive_count the number of scalars that is received by each processor
     */
#ifdef MERCURYDPM_USE_MPI
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    allGather(T& send_t, int send_count, std::vector<T>& receive_t, int receive_count)
    {
        MPI_Allgather(&send_t, send_count, Detail::toMPIType(send_t),
                                receive_t.data(), receive_count, Detail::toMPIType(receive_t[0]),communicator_);  
    }
#endif
    
    /*!
     * \brief Get the unique identifier associated with this processor.
     */
    std::size_t getProcessorID();
    
    /*!
     * \brief Get the total number of processors participating in this simulation.
     */
    std::size_t getNumberOfProcessors() const;
    
    /*!
     * \brief Get the communicator used for MPI commands.
     */
#ifdef MERCURYDPM_USE_MPI
    MPI_Comm& getComm();
#endif
    
    /*!
     * \brief Creates the MPI types telling the MPI interface how each data object looks like
     * \details NOTE: The current manner of creating MPI data types might not be compatible when computing
     * on different computers with different architecture and compilers. The "padding" which different
     * compilers and computers add to the class might be different which has a huge effect on the
     * mpi data type because it needs to be consistent over all computers. To resolve this, one
     * can easily create a consistent data type for all types required. I wish goodluck to the
     * person that needs it, because it is a lot of type work and for now I am not doing it. /MX
     */
    template<typename T>
    void createMercuryMPIType(T t, MercuryMPIType type)
    {
#ifdef MERCURYDPM_USE_MPI
        MPI_Datatype MPIType;
        MPI_Type_contiguous(sizeof(T), MPI_BYTE, &MPIType);
        MPI_Type_commit(&MPIType);
        dataTypes_.push_back(MPIType);
#endif
    }
    
    
    /*!
     * \brief Deletes the MercuryMPITypes
     */
    void deleteMercuryMPITypes()
    {
#ifdef MERCURYDPM_USE_MPI
        for(MPI_Datatype type : dataTypes_)
        {
            MPI_Type_free(&type);
        }
#endif
    }
    
    /*!
    * \brief Copy constructor is disabled, to enforce a singleton pattern
    */
    MPIContainer(const MPIContainer& orig) = delete;

private:
    
    /*!
     * \brief Constructor
     */
    MPIContainer();
    
    /*!
     * \brief The ID of the processor this class is running on
     */
    //std::size_t processorID_;
    int processorID_;
    
    /*!
     * \brief The total number of processors in the communicator
     */
    int numberOfProcessors_;

#ifdef MERCURYDPM_USE_MPI
    /*!
     * \brief List of send/receive requests that have not been resolved
     */
    std::vector<MPI_Request> pending_;
    
    /*!
     * \brief communicator that can send and recieve commands to other processors
     */
    MPI_Comm communicator_;
    /*!
     * \brief vector that contains the MercuryMPIType MPIType objects
     */ 
    std::vector<MPI_Datatype> dataTypes_;

#endif

};


#endif /* MPICONTAINER_H_ */
