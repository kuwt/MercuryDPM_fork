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

#ifndef BASEOBJECT_H
#define BASEOBJECT_H

#include <iostream>

// Foward declation of base object.
class BaseObject;

/*! 
 * \brief Operator overloading for passing the data from the BaseObject "o" into the output stream.
 */
std::ostream& operator<<(std::ostream& os, const BaseObject& o);

/*!
 * \brief Operator overloading for reading the data from an input stream into the BAseObject "o"  
 */
std::istream& operator>>(std::istream& is, BaseObject& o);

/*!
 * \class BaseObject 
 * \brief It is an abstract base class due to the purely virtual functions declared below. Even if the function is
 * purely virtual, it does not imply that it cannot have a definition. Abstract classes are useful to define a 
 * interface.
 */
class BaseObject
{
public:
    /*!
     * \brief Default constructor
     */
    BaseObject() = default;
    
    /*!
     * \brief Copy constructor, copies all the objects BaseObject contains
     */
    BaseObject(const BaseObject& p) = default;
    
    /*!
     * \brief virtual destructor
     */
    virtual ~BaseObject()  = default;
    
    /*!
     * \brief A purely virtual method with an implementation which reads the index from the stream and assigns
     * it to id_
     * \param[in] is
     */
    
    friend std::ostream& operator<<(std::ostream& os, const BaseObject& o);
    
    /*!
     * \brief 
     * \details
     */
    friend std::istream& operator>>(std::istream& is, BaseObject& o);
    
    /*!
     * \brief 
     * \param[in]
     */
    virtual void read(std::istream& is) = 0;
    
    /*!
     * \brief A purely virtual function which has an implementation which writes the name and the object id_ 
     * to the output stream.
     */
    virtual void write(std::ostream& os) const = 0;
    
    /*!
     * \brief A purely virtual function
     */
    virtual std::string getName() const = 0;
    
    /*!
     * \brief Except that it is virtual, it does the same thing as setIndex() does.
     */
    virtual void moveInHandler(unsigned int index);
    
    /*!
     * \brief Allows one to assign an index to an object in the handler/container.
     */
    void setIndex(unsigned int index);
    
    /*!
     * \brief Assigns a unique identifier to each object in the handler (container) which remains
     * constant even after the object is deleted from the container/handler.
     */
    void setId(unsigned long id);
    
    /*!
     * \brief Returns the index of the object in the handler.
     */
    unsigned int getIndex() const
    { return index_; };

    /*!
     * \brief Returns the unique identifier of any particular object.
     * \return id number of the current object
     */
    unsigned int getId() const
    { return id_; }
    
    /*!
     * \see groupId_.
     */
    void setGroupId(unsigned groupId)
    { groupId_ = groupId; }
    
    /*!
     * \see groupId_.
     */
    unsigned getGroupId() const
    { return groupId_; }

private:
    /*!
     * \brief location in BaseHandler::objects_
     */
    unsigned int index_ = 0;
    
    /*!
     * \brief unique identifier within handler (remains constant even if particle is moved)
     */
    unsigned int id_ = 0;
    
    /*!
     * \brief Identifier of a group within handler.
     * \details Useful to define properties to particle or wall clusters, etc.
     * E.g. one can apply an action to a specific group, visualise a specific group, coarse-grain a specific group, etc.
     * \see BaseHandler::nextGroupId_
     */
    unsigned int groupId_ = 0;
    
};

#endif
