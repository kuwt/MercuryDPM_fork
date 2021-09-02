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

#include <GeneralDefine.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <locale>
#include <Logger.h>

#ifndef CSVREADER_H
#define CSVREADER_H

/*!
 * \brief Enables reading of .csv files into MercuryDPM.
 *
 * \details the csvReader stores a 2D array of strings from a comma-separated value (.csv) file. In this way any data
 * type from a .csv file can be read into Mercury and converted into other formats.
 */
class csvReader
{

public:
    /*!
     * \brief Constructor
     */
    csvReader();
    
    /*!
     * \brief Reads .csv files into Mercury.
     */
    void read(const std::string& filename);
    
    /*!
     * \brief Get the Header string vector of a .csv file.
     */
    std::vector<std::string> getHeaderVector();
    
    /*!
     * \brief Set the boolean hasHeader_.
     */
    void setHeader(bool headings);
    
    /*!
     * \brief Get the 2D array with the .csv file values.
     */
    std::vector<std::vector<std::string>> getNumArray();
    
    /*!
     * \brief Get first column of a .csv file and return it as a double.
     */
    std::vector<Mdouble> getFirstColumn(Mdouble scalingFactor);
    
    /*!
     * \brief Get second column of a .csv file and return it as a double.
     */
    std::vector<Mdouble> getSecondColumn(Mdouble scalingFactor);

private:
    
    /*!
     * 2D array containing all .csv file values
     */
    std::vector<std::vector<std::string>> numArray_;
    
    /*!
     * vector of strings containing possible header values of the .csv file's first row. This vector is only filled
     * when headerFlag is set to TRUE.
     */
    std::vector<std::string> headerVector_;
    
    /*!
     * vector of doubles containing all values of the .csv file's first column.
     */
    std::vector<Mdouble> CSVFirstColumn_;
    
    /*!
     * vector of doubles containing all values of the .csv file's first column.
     */
    std::vector<Mdouble> CSVSecondColumn_;
    
    /*!
     * if FALSE the file has no headers and if TRUE the file has headers and the first row will be skipped.
     */
    bool hasHeader_ = false;
};

#endif //CSVREADER_H
