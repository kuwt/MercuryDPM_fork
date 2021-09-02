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

#include "csvReader.h"

/*!
 * \details Default constructor; sets every data member to 0 or default.
 */
csvReader::csvReader()
= default;

/*!
 * \details Read function that divides each line of the .csv file into comma delimited fields which get stored in the
 * numArray_. Further it checks for a byte order mark (BOM) and skips them if existent to avoid errors. If the
 * hasHeader_ is set to TRUE the first row of the .csv file will be skipped.
 * \param[in] filename                  Name of the .csv file.
 */
void csvReader::read(const std::string& filename)
{
    // strings to catch each line of a file and separate it into fields
    std::string line, field;
    // array of values for one line only
    std::vector<std::string> v;
    // Create ifstream object to open the file
    std::ifstream file;
    file.open(filename);
    if (file.is_open())
    {
        // get first three bytes of the ifstream object (of your file)
        int ch1 = file.get();
        int ch2 = file.get();
        int ch3 = file.get();
        // if TRUE, the file contains UTF-8 BOM (problem usually encountered on Windows OS). A BOM is a mark which
        // distinguishes different file encodings from each other. They usually start with a combination of different
        // bytes (for UTF-8-BOM it is \uFEFF, but asking for it does not work, so it is seperated into three different
        // bytes). Another workaround is changing the encoding to UTF-8 without BOM which some readers are capable of.
        if (ch1 == 0xEF && ch2 == 0xBB && ch3 == 0xBF)
        {
            // skip first three bytes if BOM exists
            file.seekg(3);
            // drop a log message to tell the user that he has a byte-order-mark
            logger(INFO, "Your file contains a byte-order-mark. The first three bytes of this file will be skipped");
        }
        else
        {
            // ensure that the ifstream object starts at the first byte of the csv file
            file.seekg(0);
        }
        while (getline(file, line))    // get next line in file
        {
            v.clear();
            std::stringstream ss(line);
            while (getline(ss, field, ','))
            {
                // break line into comma delimitted fields
                v.push_back(field);  // add each field to the 1D array
            }
            if (hasHeader_)
            {
                headerVector_ = v;  // add the 1D array to the 2D array
                hasHeader_ = false;
            }
            else
            {
                numArray_.push_back(v);  // add the 1D array to the 2D array
            }
        }
    }
    else
    {
        logger(ERROR, "File could not be opened");
    }
}

/*!
 * \details Gets the headerVector_ containing the values of the .csv file's first row.
 * \return A vector of strings containing the header values.
 */
std::vector<std::string> csvReader::getHeaderVector()
{
    return headerVector_;
}

/*!
 * \details Sets the boolean hasHeader_ which determines if the .csv file has headings.
 */
void csvReader::setHeader(bool hasHeader)
{
    hasHeader_ = hasHeader;
}

/*!
 * \details Gets the 2D array of strings which contains all the .csv file values from the read() function.
 * \return A 2D array of strings containing all the .csv file values.
 */
std::vector<std::vector<std::string>> csvReader::getNumArray()
{
    return numArray_;
}

/*!
 * \details Gets the first column of the .csv file as a vector of doubles. It extracts the first column from the 2D
 * array and converts the vector of strings into a vector of doubles.
 * \return A vector of doubles containing the values of the .csv file's first column.
 */
std::vector<Mdouble> csvReader::getFirstColumn(Mdouble scalingFactor)
{
    for (int i = 0; i < numArray_.size(); i++)
    {
        
        // converts vector<string> to vector<double>
        // if std::stod throws an exception of an invalid argument there might be an encoding issue.
        CSVFirstColumn_.push_back(std::stod(numArray_[i][0].c_str()) / scalingFactor);
    }
    return CSVFirstColumn_;
}

/*!
 * \details Gets the second column of the .csv file as a vector of doubles. It extracts the second column from the 2D
 * array and converts the vector of strings into a vector of doubles.
 * \return A vector of doubles containing the values of the .csv file's second column.
 */
std::vector<Mdouble> csvReader::getSecondColumn(Mdouble scalingFactor)
{
    for (int i = 0; i < numArray_.size(); i++)
    {
        // converts vector<string> to vector<double>
        // if std::stod throws an exception of an invalid argument there might be an encoding issue.
        CSVSecondColumn_.push_back(std::stod(numArray_[i][1].c_str()) / scalingFactor);
    }
    return CSVSecondColumn_;
}