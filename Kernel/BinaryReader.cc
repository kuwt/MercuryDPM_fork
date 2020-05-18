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



#ifndef BINARY_READER_H
#define BINARY_READER_H


#include<fstream>
#include<iostream>
#include<Logger.h>
#include<BinaryReader.h>

/*!
 * \param[in[ fileName fileName of the binary file which is to be opened
 * \deatails This is the default and only constuctor and it calls the private method openFile.
 */
BinaryReader::BinaryReader(std::string fileName)
{
    openFile(fileName);
}

/*!
 * \details This is the destructor it simple closes the file as this is the only memory that needs freeing up
 */
BinaryReader::~BinaryReader()
{
    closeFile();
}

/*!
 \details Simple function that closes the file; this is private and can only be called by the destructor
 */
void BinaryReader::closeFile()
{
    binaryFile_.close();
}

/*!
 * \param[in] fileName fileName of the binary file which is to be opened
 * \details This opens the file associated with fileName, note the file is points is stored in the private varibles binaryFile_
 */
void BinaryReader::openFile(std::string fileName)
{
    binaryFile_.open(fileName, std::ios::binary);
    if (!binaryFile_)
    {
        logger(ERROR, "BinaryReader::openFile Could not open file; file % not found", fileName);
    }
}

/*!
 * \param[in] numChar The number of characters to be read in from the binary file.
 * \return The function returns a std::string which is as long as the numChar requested.
 * \details The usages is myString = readString(6) will read in the next 6 characters from the binaryfile as a string.
 */
std::string BinaryReader::readString(unsigned int numChar)
{
    char tempRead[numChar];
    binaryFile_.read(tempRead, numChar);
    std::string returnString(tempRead);
    return returnString;
    
}

/*!
 * \param[in] size This is the size in bytes of the double stored in the file
 * \returns a double that is the double read from the file
 * \details Usage myDouble = readDouble(4) would read a 8*4 = 32 bit double from the file; or myDouble = readDouble(8) would read a 64 bit double from the file.
 */
double BinaryReader::readDouble(unsigned int size)
{
    char tempRead[size];
    binaryFile_.read(tempRead, 4);
    double* returnDoublePtr = reinterpret_cast<double*>(tempRead);
    return (*returnDoublePtr);
}

/*!
 * \param[in] size This is the size in bytes of the float stored in the file
 * \returns a double (not a float) that is the upcast of the float that is read in.
 * \details Usage myDouble = readFloat(4) would read 32 bit float from file covert it to a double and then return this double. \see BinaryReader::readDouble
 * Note, this is different to readDouble as the binary represation of a double and a float is different.
 */
double BinaryReader::readFloat(unsigned int size)
{
    char tempRead[size];
    binaryFile_.read(tempRead, 4);
    float* returnFloatPtr = reinterpret_cast<float*>(tempRead);
    return static_cast<double>(*returnFloatPtr);
}

/*!
 * \param[in] size This is the size in bytes of the unsigned int that is stored in the file.
 * \returns a unsigned int this is the unisgned int that is read in
 * \details Usage myUnsignedInt = readUnsignedInt(4) woudl read a 32 bit (4*8) unsigned int from the file.
 */
unsigned int BinaryReader::readUnsignedInt(unsigned int size)
{
    char tempRead[size];
    binaryFile_.read(tempRead, size);
    unsigned int* returnUnsignedInt = reinterpret_cast<unsigned int*>(tempRead);
    return (*returnUnsignedInt);
}

/*!
 * \param[in] size number of bytpes to be ignored
 * \details Reads and disgards the next size characters size*8 bytes from the binary file
 */
void BinaryReader::ignoreChar(unsigned int size)
{
    char tempRead[size];
    binaryFile_.read(tempRead, size);
}

#endif
