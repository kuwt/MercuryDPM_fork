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
 * \file BaseHandler.h
 *
 * Class BaseHandler
 */
#ifndef BASEHANDLER_H
#define BASEHANDLER_H

#include <vector>
#include <type_traits>
#include "Math/Helpers.h"
#include "Logger.h"

class DPMBase;

/*!
 * \class BaseHandler
 * \brief Container to store the pointers to all objects that one creates in a simulation. 
 * \details The BaseHandler allows one to create a container to store all pointer objects of a templated type T 
 * It is implemented by a (protected) vector of pointers to objects of type T. Once the container is created, the BaseHandler
 * also provides the provision to manipulate the pointers i.e. by accessing, adding, deleting and few more operations by using its
 * member methods. 
 */
template<typename T>
class BaseHandler
{
public:
    /*!
     * \brief Default BaseHandler constructor, it creates an empty BaseHandler and assigns DPMBase_ to a
     * null pointer.
     */
    BaseHandler();
    
    /*!
     * \brief Constructor that copies the objects of the given handler into itself and sets other variables to 0/nullptr.
     */
    BaseHandler(const BaseHandler<T>& BH);
    
    /*!
     * \brief Destructor, it destructs the BaseHandler and all Object it contains.
     */
    virtual ~BaseHandler();
    
    /*!
     * \brief Function that copies the contents (vector of pointers, maxObject_, nextId_, DPMBase_) from one 
     * handler (container) to the other.
     */
    void copyContentsFromOtherHandler(const BaseHandler<T>& BH);
    
    /*!
     * \brief Creates a copy of a Object and adds it to the BaseHandler.
     */
    template<typename U>
    typename std::enable_if<!std::is_pointer<U>::value, U*>::type
    copyAndAddObject(const U& object);
    
    /*!
     * \brief Creates a copy of a Object and adds it to the BaseHandler.
     */
    template<typename U>
    typename std::enable_if<std::is_pointer<U>::value, U>::type
    copyAndAddObject(const U object);
    
    /*!
     * \brief Creates a copy of a Object and adds it to the BaseHandler.
     * This is one locally for inserting mpi particles, they avoid the global
     * check if the particle can actually be inserted, because the mpi domain already
     * knows that is the case
     */
    template<typename U>
    typename std::enable_if<!std::is_pointer<U>::value, U*>::type
    copyAndAddGhostObject(const U& object);
    
    /*!
     * \brief Creates a copy of a Object and adds it to the BaseHandler.
     * This is one locally for inserting mpi particles, they avoid the global
     * check if the particle can actually be inserted, because the mpi domain already
     * knows that is the case
     */
    template<typename U>
    typename std::enable_if<std::is_pointer<U>::value, U>::type
    copyAndAddGhostObject(const U object);
    
    /*!
     * \brief Adds an existing object to the BaseHandler without 
     * changing the id of the object
     */
    virtual void addExistingObject(T* O);
    
    /*!
     * \brief Adds a new Object to the BaseHandler.
     */
    virtual void addObject(T* object);
    
    /*!
     * \brief Adds a new Object to the BaseHandler.
     * called by the to avoid increasing the id
     */
    virtual void addGhostObject(T* O);
    
    void removeIf(const std::function<bool(T*)> cond);
    
    /*!
     * \brief Removes an Object from the BaseHandler.
     */
    virtual void removeObject(unsigned const int index);
    
    /*!
     * \brief Removes the last Object from the BaseHandler.
     */
    void removeLastObject();
    
    /*!
     * \brief Empties the whole BaseHandler by removing all Objects and setting all other variables to 0.
     */
    virtual void clear();
    
    /*!
     * \brief Reads Object into the BaseHandler from restart data.
     * \param[in] is The input stream from which the information is read.
     */
    virtual void readAndAddObject(std::istream& is) = 0;
    
    /*!
     *  \brief Reads all objects from restart data.
     */
    void read(std::istream& is);
    
    /*!
     * \brief Gets a pointer to the Object at the specified index in the BaseHandler. 
     */
    T* getObjectById(const unsigned int id);
    
    /*!
     * \brief Gets a vector of pointers to the objects with the specific id.
     */
    std::vector<T*> getObjectsById(const unsigned int id);
    
    /*!
     * \brief Gets a pointer to the Object at the specified index in the BaseHandler.   
     */
    T* getObject(const unsigned int id);
    
    /*!
     * \brief Gets a constant pointer to the Object at the specified index in the BaseHandler.
     */
    const T* getObject(const unsigned int id) const;
    
    /*!
     * \brief Gets a pointer to the last Object in this BaseHandler.
     */
    T* getLastObject();
    
    /*!
     * \brief Gets a constant pointer to the last Object in this BaseHandler.
     */
    const T* getLastObject() const;
    
    /*!
     * \brief Gets the number of real Object in this BaseHandler. (i.e. no mpi or periodic particles)
     */
    virtual unsigned int getNumberOfObjects() const;
    
    /*!
     * \brief Gets the size of the particleHandler (including mpi and periodic particles)
     */
    unsigned int getSize() const;
    
    /*!
     * \brief Gets the storage capacity of this BaseHandler.
     */
    unsigned int getStorageCapacity() const;
    
    /*!
     * \brief Sets the storage capacity of this BaseHandler.
     */
    void setStorageCapacity(const unsigned int N);
    
    /*!
     * \brief Resizes the container to contain N elements
     */
    void resize(const unsigned int N, const T& obj);
    
    /*!
     * \brief Gets the begin of the const_iterator over all Object in this BaseHandler.
     */
    const typename std::vector<T*>::const_iterator begin() const;
    
    /*!
     * \brief Gets the begin of the iterator over all BaseBoundary in this BaseHandler.
     */
    const typename std::vector<T*>::iterator begin();
    
    /*!
     * \brief Gets the end of the const_iterator over all BaseBoundary in this BaseHandler.
     */
    const typename std::vector<T*>::const_iterator end() const;
    
    /*!
     * \brief Gets the end of the iterator over all BaseBoundary in this BaseHandler.
     */
    const typename std::vector<T*>::iterator end();
    
    /*!
     * \brief Sets the problem that is solved using this handler.
     */
    /// \todo MX: Bad practice to have a function parameter with the exact name of the class
    void setDPMBase(DPMBase* DPMBase);
    
    /** This function sets the id and ensures that nextId is a bigger value than id. 
     * \todo we should use this function only to set the id of particles, not BaseObject::setId; however, to block BaseObject::setId, I need to make this function a friend of BaseObject, and I don't know how to do that.
     */
    void setId(T* object, unsigned int id)
    {
        object->setId(id);
        if (nextId_ <= id)
        {
            nextId_ = id + 1;
        }
    }
    
    /*
     * \brief This function updates the iD counter by one
     * \details Function used in parallel to keep the Id's unique over all processors
     */
    void increaseId()
    {
        nextId_++;
    }
    
    unsigned int getNextId()
    {
        return nextId_;
    }
    
    void setNextId(unsigned int id)
    {
        nextId_ = id;
    }
    
    /*! 
     * \brief Gets the problem that is solved using this handler.
     */
    DPMBase* getDPMBase();
    
    /*! 
     * \brief Gets the problem that is solved using this handler and does not change the class.
     */
    DPMBase* getDPMBase() const;
    
    /*!\brief Gets the name of this handler.
     * \return A string that contains the name of the handler.
     */
    virtual std::string getName() const = 0;
    
    /*!
     * \brief now empty function for writing VTK files.
     * \deprecated Now the VTK-writers have their own classes, and are called from DPMBase.
     */
    [[deprecated]]
    virtual void writeVTK() const
    {};
    
    /*!
     * \brief Should be called each time you assign a groupId.
     * Returns the value of nextGroupId_ and increases nextGroupId_ by one.
     */
    unsigned getNextGroupId() { return nextGroupId_++; }
    
protected:
    /*!
     * \brief The actual list of Object pointers
     * 
     * The list of Object pointers. This handler is responsible for the memory-deallocation
     * of these objects.
     */
    std::vector<T*> objects_;

private:
    /*!
     * \brief An integer to keep track of the largest number of objects ever stored in this BaseHandler
     */
    unsigned int maxObjects_;
    
    /*!
     * \brief identifier for next object created
     */
    unsigned int nextId_;
    
    /*!
     * \brief value of the next BaseObject::groupId_.
     * Value increased by one each time a groupId is assigned.
     * Default group is 0
     */
    unsigned nextGroupId_ = 1;

    /*!
     * \brief A pointer back to the DPMBase class.
     * 
     * Please note that this pointer back to the DPMBase class is a "shared" pointer
     * and should not be deallocated by this class.
     */
    DPMBase* DPMBase_;
};

/*!
 * 
 */
template<typename T>
BaseHandler<T>::BaseHandler()
{
    DPMBase_ = nullptr;
    clear();
    logger(DEBUG, "BaseHandler<T>::BaseHandler() finished");
}

/*!
 * \param[in] BH A reference to the BaseHandler that has to be copied.
 * \details This is not a copy constructor! It only copies the vector objects_
 *          from the given handler, and sets all other variables to 0/nullptr.
 * \todo Should max objects be set to the number of objects after this constructor?
 *       Maybe in copyContentsFromOtherHandler?
 */
template<typename T>
BaseHandler<T>::BaseHandler(const BaseHandler<T>& BH)
{
    DPMBase_ = nullptr;
    clear();
    copyContentsFromOtherHandler(BH);
    logger(DEBUG, "BaseHandler<T>::BaseHandler(const BaseHandler &BH) finished");
}

template<typename T>
BaseHandler<T>::~BaseHandler()
{
    clear();
    logger(DEBUG, "BaseHandler<T>::~BaseHandler() finished");
}

///\param[in] BH A reference to the BaseHandler of which the objects have to be copied.
template<typename T>
void BaseHandler<T>::copyContentsFromOtherHandler(const BaseHandler<T>& BH)
{
    for (const T* const obj : BH.objects_)
    {
        addObject(obj->copy());
    }
}

///\param[in] object A reference to the BaseHandler of which the objects have to be copied.
template<typename T>
template<typename U>
typename std::enable_if<!std::is_pointer<U>::value, U*>::type
BaseHandler<T>::copyAndAddObject(const U& object)
{
    U* oCopy = object.copy();
    addObject(oCopy);
    return oCopy;
}

///\param[in] object A reference to the Object that has to be copied.
template<typename T>
template<typename U>
typename std::enable_if<std::is_pointer<U>::value, U>::type
BaseHandler<T>::copyAndAddObject(const U object)
{
    return copyAndAddObject(*object);
}

///\param[in] object A reference to the BaseHandler of which the objects have to be copied.
template<class T>
template<class U>
typename std::enable_if<!std::is_pointer<U>::value, U*>::type
BaseHandler<T>::copyAndAddGhostObject(const U& object)
{
    U* oCopy = object.copy();
    addGhostObject(oCopy);
    return oCopy;
}

///\param[in] object A reference to the Object that has to be copied.
template<class T>
template<class U>
typename std::enable_if<std::is_pointer<U>::value, U>::type
BaseHandler<T>::copyAndAddGhostObject(const U object)
{
    return copyAndAddGhostObject(*object);
}

///\param[in] object A point to an existing object which already has an id given
template<class T>
void BaseHandler<T>::addExistingObject(T* O)
{
    objects_.push_back(O);
    //Set the index of the particle
    getLastObject()->setIndex(getSize() - 1);
    //Adjust the nextId_ value
    if (O->getId() + 1 > nextId_)
    {
        nextId_ = O->getId() + 1;
    }
}

///\param[in] object A pointer to the object that must be added.
template<class T>
void BaseHandler<T>::addObject(T* object)
{
    objects_.push_back(object);
    //Set the index of the particle
    getLastObject()->setIndex(getSize() - 1);
    //set the non changing particle identifier
    getLastObject()->setId(nextId_);
    //Update Id for next particle
    nextId_++;
}

/// \todo mx: type the stuff here: keeps the id unique key
template<class T>
void BaseHandler<T>::addGhostObject(T* O)
{
    objects_.push_back(O);
    //Set the index of the particle
    getLastObject()->setIndex(getSize() - 1);
}

//\todo should this function ever be used?
template<class T>
void BaseHandler<T>::removeIf(const std::function<bool(T*)> cond)
{
    for (int i = 0; i < objects_.size(); ++i) {
        if (cond(objects_[i]))
        {
            //objects_(i)->actionsOnErase();
            removeObject(i);
            --i;
        }
    }
}

/*!
 * This methods removes an object. This methods invalidates ANY iterators to
 * objects in this container. This method may shuffle the order of objects in
 * this container.
 * \param[in] index An unsigned integer that gives the id of the Object that has to be removed.
 */
template<typename T>
void BaseHandler<T>::removeObject(const unsigned int index)
{
    logger.assert(index < getSize(),
                  "In: void %::removeObject(const unsigned int index) const, "
                  "no object exists with index %, number of objects is %",
                  getName(), index, getSize());
    
    //Removing a particle within the list is not efficient.
    //Swap with last element and then perform the deletion
    const unsigned int lastIndex = objects_.size() - 1;
    
    T* const objectToDelete = objects_[index];
    
    //No swapping required if it is the last object
    if (index != lastIndex)
    {
        
        T* const objectToMove = objects_[lastIndex];
        
        objects_[index] = objectToMove; //place it back
        objects_[lastIndex] = objectToDelete; //Just to make sure.
        
        //and notify it of the change.
        objects_[index]->moveInHandler(index);
        //Even though we are going to delete this particle,
        //we still need to keep it consistent.
        objects_[lastIndex]->moveInHandler(lastIndex);
    }
    
    
    //And clear it from the backing container.
    objects_.pop_back();
    //And _NOW_ we delete it.
    
    delete objectToDelete;
    
}

template<typename T>
void BaseHandler<T>::removeLastObject()
{
    if (getSize() == 0)
    {
        logger(WARN, "In: void %::removeLastObject, no Object exists in this BaseHandler.", getName());
        return;
    }
    T* const object = objects_.back();
    //Remove the (now double) reference to that last Object
    objects_.pop_back();
    //Physically removes Object
    delete object;
}

///Delete all objects stored in objects_ and set the maximum number of objects that
///have been in this container to 0, and set the Id of the next object that will be added to 0.
template<typename T>
void BaseHandler<T>::clear()
{
    
    for (T* const obj : objects_)
    {
        delete obj;
    }
    objects_.clear();
    
    nextId_ = 0;
    maxObjects_ = 0;
}

/// \param[in] is The input stream from which the information is read.
template<typename T>
void BaseHandler<T>::read(std::istream& is)
{
    clear();
    unsigned int N;
    std::string dummy;
    is >> dummy;
    std::stringstream line;
    helpers::getLineFromStringStream(is, line);
    line >> N;
    logger(VERBOSE, "In %::read(is): reading in % objects.", getName(), N);
    setStorageCapacity(N);
    for (unsigned int i = 0; i < N; i++)
    {
        readAndAddObject(is);
    }
}

/// \param[in] id The id of the requested Object.
/// \return A pointer to the Object with the correct Id.   
///Gets an object with the identity id. Please note that the object with this identity
///does not have to be at place id in the vector of Object objects_.
template<typename T>
T* BaseHandler<T>::getObjectById(const unsigned int id)
{
    // Usually, the id and the index into the backing storage matches
    // So check this position first!
    // dducks: Can't we guarantee more? That should speed up searches. 
    if (id < objects_.size() && objects_[id]->getId() == id)
    {
        return objects_[id]; //There is a hit, return early
    }
    
    for (T* obj : objects_) //Search for the correct id, since it wasn't where
    {                        // we expected it. Just use a linear search..
        if (obj->getId() == id) //Found it, so return!
            return obj;
    }
#ifndef MERCURY_USE_MPI
    logger(ERROR, "[BaseHandler::getObjectById()] in Object* %: Object with ID % could not be found.", getName(), id);
#endif
    return nullptr;
}

/*!
 * \details Gets all the objects with the id from the handler. This is especially useful in parallel
 * where ghost particles of a real particle have the same id
 * \param[in] id The id of the requested objects.
 * \return A list of pointers towards objects with the specified id
 */
template<typename T>
std::vector<T*> BaseHandler<T>::getObjectsById(const unsigned int id)
{
    std::vector<T*> list;
    for (T* obj : objects_)
    {
        if (obj->getId() == id)
        {
            list.push_back(obj);
        }
    }
#ifndef MERCURY_USE_MPI
    logger(ERROR, "[BaseHandler::getObjectById()] in Object* %: Object with ID % could not be found.", getName(), id);
#endif
    return list;
    
}

///\param[in] index the index of the requested Object.
///\return A pointer to the requested Object.  
template<typename T>
T* BaseHandler<T>::getObject(const unsigned int index)
{
    logger.assert(index < getSize(),
                  "[%::getObject()] Object couldn't be found because index (%) is higher than number of objects (%).",
                  getName(), index, getSize());
    return objects_[index];
}

/// \param[in] index the index of the requested Object.
/// \return A constant pointer to the requested Object.
template<typename T>
const T* BaseHandler<T>::getObject(const unsigned int index) const
{
    logger.assert(index < getSize(),
                  "[%::getObject() const] Object couldn't be found because index (%) is higher than number of objects (%).",
                  getName(), index, getSize());
    return objects_[index];
}

///\return A pointer to the last Object in the BaseHandler.
template<typename T>
T* BaseHandler<T>::getLastObject()
{
    return objects_.back();
}

///\return A constant pointer to the last Object in the BaseHandler.
template<typename T>
const T* BaseHandler<T>::getLastObject() const
{
    return objects_.back();
}

///\return The number of Objects in this BaseHandler.
template<typename T>
unsigned int BaseHandler<T>::getNumberOfObjects() const
{
    return objects_.size();
}

///\return The number of items in this BaseHandler
template<class T>
unsigned int BaseHandler<T>::getSize() const
{
    return objects_.size();
}

///\return The storage capacity of this BaseHandler.
template<typename T>
unsigned int BaseHandler<T>::getStorageCapacity() const
{
    return objects_.capacity();
}

///\param[in] N The storage capacity the BaseHandler will have 
template<typename T>
void BaseHandler<T>::setStorageCapacity(const unsigned int N)
{
    objects_.reserve(N);
}

///If the current size is greater than count, the container is reduced to its first count elements.
///If the current size is less than count, additional elements are appended and initialized with copies of obj.
///\param[in] N Container resized to contain N elements.
///\param[in] obj additional elements are appended and initialized with copies of obj.
template<class T>
void BaseHandler<T>::resize(const unsigned int N, const T& obj)
{
    //objects_.resize(N,obj); //doesn't work because the handler stores pointers only (data needs to be allocated);
    while (getSize() < N)
        copyAndAddObject(obj);
    while (getSize() > N)
        removeLastObject();
}

///\return A const_iterator pointing to the first Object.
template<typename T>
const typename std::vector<T*>::const_iterator BaseHandler<T>::begin() const
{
    return objects_.begin();
}

///\return A iterator pointing to the first Object.
template<typename T>
const typename std::vector<T*>::iterator BaseHandler<T>::begin()
{
    return objects_.begin();
}

///\return A const_iterator pointing to the last BaseBoundary.
template<typename T>
const typename std::vector<T*>::const_iterator BaseHandler<T>::end() const
{
    return objects_.end();
}

///\return An iterator pointing to the last BaseBoundary.
template<typename T>
const typename std::vector<T*>::iterator BaseHandler<T>::end()
{
    return objects_.end();
}

///\param[in] DPMBase A pointer to a DPMBase, which is the superclass for all problem descriptions.
template<typename T>
void BaseHandler<T>::setDPMBase(DPMBase* DPMBase)
{
    DPMBase_ = DPMBase;
}

///\return A pointer to the DPMBase (problem descriptor) that is using this handler.
template<typename T>
DPMBase* BaseHandler<T>::getDPMBase()
{
    return DPMBase_;
}

///\return A pointer to the DPMBase (problem descriptor) that is using this handler.
template<typename T>
DPMBase* BaseHandler<T>::getDPMBase() const
{
    return DPMBase_;
}

#endif

