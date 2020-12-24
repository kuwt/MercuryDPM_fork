//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://MercuryDPM.org/Team>.
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

#ifndef BinaryReader_H
#define BinaryReader_H
#include<fstream>
#include<iostream>
#include<string>

/*!
 * \brief This gives functionality to read information from binary formats like STL etc.
 * This class is complete stand-alone and is tested with one any reference to other MecuryDPM code except Vections and Logger.
 */
class BinaryReader
{
public:
    
    /*!
     * \brief Default constuction, requires to users to prove the name of the file that will be opened.
     */
    explicit BinaryReader(std::string);
    
    /*!
     * \brief Destructor, simple closes the file
     */
    ~BinaryReader();
    
    /*!
     * \brief reads the next so many Characters (bytes) as a std::string
     */
    std::string readString(unsigned int numChar);
    
    /*!
     * \brief read the next so many bytes as a double
     */
    double readDouble(unsigned int size);
    
    /*!
     * \brief read the next so many bytes as a unsined int
     */
    unsigned int readUnsignedInt(unsigned int size);
    
    /*!
     * \brief read the next so many bytes as a double (not in this case they were saves as a float orgainlly)
     */
    double readFloat(unsigned int size);
    
    /*!
     * \brief read and ignore the next number of characters
     */
    void ignoreChar(unsigned int size);

private:
    
    /*!
     * \brief opens the file with fileName
     */
    void openFile(std::string fileName);
    
    /*!
     * \brief close the file with fileName
     */
    void closeFile();
    
    /// The pointer for the binary file.
    std::ifstream binaryFile_;
    
};
#endif