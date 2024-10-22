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


#include "File.h"

#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "Logger.h"

/*!
 * \brief Pads the number
 * This function tries to pad the number to 4 digits, which is used when
 * you create multiple files with padded numbers. Any numbers larger than
 * 4 digits return unmodified.
 * \param value The value to modify
 * \returns A padded string
 */
std::string to_string_padded(unsigned int value)
{
    std::ostringstream out;
    out << std::setw(4) << std::setfill('0') << value;
    return out.str();
}

/*!
 * \param[in,out] os output stream to which the fileType is written
 * \param[in] fileType the fileType that has to be written to the output stream
 * \return the output stream "os" that is returned after adding the fileType string
 */
std::ostream& operator<<(std::ostream& os, FileType fileType)
{
    if (fileType == FileType::NO_FILE)
        os << "NO_FILE";
    else if (fileType == FileType::ONE_FILE)
        os << "ONE_FILE";
    else if (fileType == FileType::MULTIPLE_FILES)
        os << "MULTIPLE_FILES";
    else if (fileType == FileType::MULTIPLE_FILES_PADDED)
        os << "MULTIPLE_FILES_PADDED";
    else
    {
        logger(ERROR, "FileType not recognized");
    }
    return os;
}

/*!
 * \param[in,out] is The input stream from which the fileType is read
 * \param[in] fileType The fileType that has to be read from the input stream
 * \return the input stream "is" (that is returned after the fileType string is read out)
 */
std::istream& operator>>(std::istream& is, FileType& fileType)
{
    std::string fileTypeString;
    is >> fileTypeString;
    if (!fileTypeString.compare("NO_FILE"))
        fileType = FileType::NO_FILE;
    else if (!fileTypeString.compare("ONE_FILE"))
        fileType = FileType::ONE_FILE;
    else if (!fileTypeString.compare("MULTIPLE_FILES"))
        fileType = FileType::MULTIPLE_FILES;
    else if (!fileTypeString.compare("MULTIPLE_FILES_PADDED"))
        fileType = FileType::MULTIPLE_FILES_PADDED;
    else
    {
        logger(ERROR, "operator>>: FileType % not recognized", fileTypeString);
    }
    return is;
}

/*!
 * \details A File constructor which initialises FILE::saveCount_=0, sets the default filename as FILE::name_ = "", sets the 
 * default FileType to be as FileType::ONE_FILE and few other variables as listed below.   
 */
File::File()
{
    //sets the default for the number of time steps to be skipped
    //in between each saved "snapshot" of the system to zero
    //(i.e. records every time step by default)
    saveCount_ = 0;
    
    // file name has to be set by the user
    name_ = "out";
    
    // output into a single file by default
    fileType_ = FileType::ONE_FILE;
    
    // counter of currently open file set to 0 by default
    counter_ = 0;
    
    //stores the time step of the last write/read operation; NEVER by default
    lastSavedTimeStep_ = NEVER;
    
    //sets the default openMode to "out"
    //i.e. files will by default be written to, not read from.
    openMode_ = std::fstream::out;
    
    //sets the default logarithmicSaveCount to zero
    logarithmicSaveCountBase_ = 0;
}

/*!
 * \details Copy everything but the fstream object (which cannot be copied).
 * \param[in] f the file object that is to be copied.
 */
File::File(const File& f)
{
    saveCount_ = f.saveCount_;
    name_ = f.name_;
    fileType_ = f.fileType_;
    counter_ = f.counter_;
    lastSavedTimeStep_ = f.lastSavedTimeStep_;
    openMode_ = f.openMode_;
    logarithmicSaveCountBase_ = f.logarithmicSaveCountBase_;
}

/*!
 * \details Destructor
 */
File::~File()
= default;

/*!
 * \details Returns fstream (file stream) of any file for input and output tasks.
 * \return std::fstream&
 */
std::fstream& File::getFstream()
{
    return fstream_;
}

//void File::setFstream(const std::fstream& file)
//{
//    this->fstream_ = file;
//}
/*!
 * \return Returns a constant reference of type const std::string&
 */
const std::string& File::getName() const
{
    return name_;
}

const std::string File::getFullName() const
{
    return getFullName(getCounter() - 1);
}

/*!
 * \details In case of FileType:fileType_== multiple files or multiple files padded, multiple files are generated and are named as problem.data.0, problem.data.1 or problem.data.0000, problem.data.0001 
 * \return Returns a constant of type std::string
 */
const std::string File::getFullName(unsigned counter) const
{
    //get the full file name
    std::stringstream lastName("");
    lastName << name_;
    if (getFileType() == FileType::MULTIPLE_FILES)
    {
        lastName << "." << counter;
    }
    else if (getFileType() == FileType::MULTIPLE_FILES_PADDED)
    {
        lastName << "." << to_string_padded(counter);
    }
    return lastName.str();
}

/*!
 * \param[in] name (Takes in the to be name of the File)
 */
void File::setName(const std::string& name)
{
    logger.assert_always(!getName().empty(), "Error: Name cannot be empty");
    this->name_ = name;
}

/*!
 * \return Returns the FileType (File::fileType_)
 */
FileType File::getFileType() const
{
    return fileType_;
}

/*!
 * \param[in] fileType
 */
void File::setFileType(FileType fileType)
{
    fileType_ = fileType;
}

/*!
 * \return unsigned int counter_
 */
unsigned int File::getCounter() const
{
    return counter_;
}

/*!
 * \param[in] counter
 */
void File::setCounter(unsigned int counter)
{
    counter_ = counter;
}

/*!
 * \return std::fstream::openmode
 */
std::fstream::openmode File::getOpenMode() const
{
    return openMode_;
}

/*!
 * \param[in] openmode
 */
void File::setOpenMode(std::fstream::openmode openMode)
{
    openMode_ = openMode;
}

/*!
 * \return unsigned int saveCount_
 */
unsigned int File::getSaveCount() const
{
    return saveCount_;
}

/*!
 * \details File::setSaveCount assigns the number of time steps to be skipped before the data is written to an existing file or a new file.
 * \param[in] saveCount
 */
void File::writeFirstAndLastTimeStep()
{
    saveCount_ = NEVER;
}

/*!
 * \details File::setSaveCount assigns the number of time steps to be skipped before the data is written to an existing file or a new file.
 * \param[in] saveCount 
 */
void File::setSaveCount(unsigned int saveCount)
{
    saveCount_ = saveCount;
}


/*!
 * \details File::setlogarithmicSaveCount assigns the base of the saveCount on a logarithmic time scale
 * \param[in] logarithmicSaveCountBase
 */
void File::setlogarithmicSaveCount(const Mdouble logarithmicSaveCountBase)
{
    logger.assert_always(logarithmicSaveCountBase > 1, "logarithmicSaveCountBase should always be larger than 1");
    logarithmicSaveCountBase_ = logarithmicSaveCountBase;
}

/*!
 * \details Returns the time step at which the next write or read operation has to happen
 * \return unsigned int nextSaveTimeStep_
 */
unsigned int File::getLastSavedTimeStep() const
{
    return lastSavedTimeStep_;
}

/*!
 * \details Allows one to set the time step at which the next write or read operation has to happen
 * \param[in] lastSavedTimeStep
 */
void File::setLastSavedTimeStep(unsigned int lastSavedTimeStep)
{
    lastSavedTimeStep_ = lastSavedTimeStep;
}

/*!
 * \details 
 * \param[in] ntimeSteps
 * \return True or False (a bool)
 */
bool File::saveCurrentTimeStep(unsigned int ntimeSteps) {
    return getFileType() != FileType::NO_FILE && saveCurrentTimeStepNoFileTypeCheck(ntimeSteps);
}


bool File::saveCurrentTimeStepNoFileTypeCheck(unsigned int ntimeSteps)
{
    /* check:
     * - if this time step should be written
     * - if the file type is not NO_FILE
     * - if file can be opened
     * in that case, change lastSavedTimeStep and return true;
     */
    if ((lastSavedTimeStep_ == NEVER || ntimeSteps >= lastSavedTimeStep_ + saveCount_))
    {
        //note: do not do the following at t = 0 because this makes no sense for a logarithm
        if (logarithmicSaveCountBase_ > 1 && ntimeSteps > 0 &&
            saveCount_ < ceil((logarithmicSaveCountBase_ - 1) * lastSavedTimeStep_))
        {
            /*calculate the new saveCount base on the user input logarithmicSaveCountBase,
             *and multiply by the actual number of time steps
             */
            saveCount_ = ceil((logarithmicSaveCountBase_ - 1) * lastSavedTimeStep_);
        }
        return true;
    } else {
        return false;
    }
}

/*!
 * \details Returns a bool to check if the file is open or closed. It also increments the nextSavedTimeStep with the saveCount
 * \return bool (True or False)
 * \bug Deepak checked by using fstream_.fail() instead of !fstrea_.is_open(), however this breaks selftests, Thomas will look at this
 */
bool File::open()
{
    //close old file if multi-file output
    if (fileType_ == FileType::MULTIPLE_FILES || fileType_ == FileType::MULTIPLE_FILES_PADDED)
    {
        fstream_.close();
    }
    
    counter_++;
    
    if (!fstream_.is_open())
    {
        fstream_.open(getFullName().c_str(), openMode_);
        if (!fstream_.is_open())
        {
            return false;
        }
    }
    
    return true;
}

/*!
 * \param[in] openMode
 */
bool File::open(std::fstream::openmode openMode)
{
    setOpenMode(openMode);
    return open();
}

/*!
 * \param[in] openMode
 */
bool File::openWrite(unsigned nTimeSteps)
{
    setLastSavedTimeStep(nTimeSteps);
    if (getFileType() == FileType::ONE_FILE && getCounter() != 0)
    {
        setOpenMode(std::fstream::out | std::fstream::app);
    }
    else
    {
        setOpenMode(std::fstream::out);
    }
    return open();
}

/*!
 * \param[in] openMode
 */
bool File::openWriteNoAppend(unsigned nTimeSteps)
{
    setLastSavedTimeStep(nTimeSteps);
    return open(std::fstream::out);
}

/*!
 *
 */
void File::close()
{
    fstream_.close();
}

/*!
 * \details Read function, which accepts an input stream object as input and assigns the member variables i.e. name_, fileType_,
 * saveCount_, counter_ and lastSavedTimeStep_
 * \param[in,out] is
 */
void File::read(std::istream& is)
{
    std::string dummy;
    is >> dummy;
    if (!dummy.compare("name"))
        is >> name_ >> dummy;
    is >> fileType_;
    is >> dummy >> saveCount_;
    is >> dummy >> counter_;
    if (counter_ != 0)
    {
        is >> dummy >> lastSavedTimeStep_;
    }
    else
    {
        lastSavedTimeStep_ = NEVER;
    }
    //if (dummy != "lastSavedTimeStep") lastSavedTimeStep_=SAVE;
}

/*!
 * \details BaseParticle print function, which accepts an output stream object as input 
 *  and writes the info to the std::ostream
 * \param[in,out] os
 */
void File::write(std::ostream& os) const
{
    //only write name if it differs from the default name
    if (getFullName().compare(name_))
        os << "name " << name_ << " ";
    os << "fileType " << fileType_;
    os << " saveCount " << saveCount_;
    os << " counter " << counter_;
    if (counter_ != 0)
    {
        os << " lastSavedTimeStep " << lastSavedTimeStep_;
    }
    ///\todo TW: openMode_ is not saved, maybe it should not even be stored but set every time you open a file
}

/*!
 * \param[in,out] os 
 * \param[in] o
 * \return std::ostream& os
 */
std::ostream& operator<<(std::ostream& os, const File& o)
{
    o.write(os);
    return os;
}

/*!
 * \param[in,out] is
 * \param[in] o
 * \return std::istream&
 */
std::istream& operator>>(std::istream& is, File& o)
{
    o.read(is);
    return (is);
}
