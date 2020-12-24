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


#include "BaseObject.h"
#include "Logger.h"
#include "Math/Helpers.h"

/*!
 * \param[in,out] os
 * \param[in] o
 * \return std::ostream &
 */
std::ostream& operator<<(std::ostream& os, const BaseObject& o)
{
    o.write(os);
    return os;
}

/*!
 * \param[in] o
 * \param[in,out] is
 * \return std::istream&
 */
std::istream& operator>>(std::istream& is, BaseObject& o)
{
    o.read(is);
    return (is);
}

/*!
 * \param[in] index 
 */
void BaseObject::moveInHandler(const unsigned int index)
{
    index_ = index;
}

/*!
 * \param[in] index index number of the current object
 */
void BaseObject::setIndex(const unsigned int index)
{
    index_ = index;
}

/*!
 * \param[in] id id number of the current object
 */
void BaseObject::setId(unsigned long id)
{
    id_ = id;
    ///\todo TW: here we should update BaseHandler::nextId_
}

/*!
 * \param[in] is stream object from which data is read
 */
void BaseObject::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> id_;
    helpers::readOptionalVariable(is, "groupID", groupId_);
}

/*!
 * \param[in] os stream object to which data is written
 */
void BaseObject::write(std::ostream& os) const
{
    os << getName();
    os << " id " << id_;
    if (groupId_) os << " groupID " << groupId_;
}

