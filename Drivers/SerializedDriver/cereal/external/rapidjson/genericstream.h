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

// Generic*Stream code from https://code.google.com/p/rapidjson/issues/detail?id=20
#ifndef RAPIDJSON_GENERICSTREAM_H_
#define RAPIDJSON_GENERICSTREAM_H_

#include "rapidjson.h"
#include <iostream>

namespace rapidjson {

  //! Wrapper of std::istream for input.
  class GenericReadStream {
    public:
      typedef char Ch;    //!< Character type (byte).

      //! Constructor.
      /*!
        \param is Input stream.
        */
      GenericReadStream(std::istream & is) : is_(&is) {
      }


      Ch Peek() const {
        if(is_->eof()) return '\0';
        return static_cast<char>(is_->peek());
      }

      Ch Take() {
        if(is_->eof()) return '\0';
        return static_cast<char>(is_->get());
      }

      size_t Tell() const {
        return (int)is_->tellg();
      }

      // Not implemented
      void Put(Ch)       { RAPIDJSON_ASSERT(false); }
      void Flush()       { RAPIDJSON_ASSERT(false); }
      Ch* PutBegin()     { RAPIDJSON_ASSERT(false); return 0; }
      size_t PutEnd(Ch*) { RAPIDJSON_ASSERT(false); return 0; }

      std::istream * is_;
  };


  //! Wrapper of std::ostream for output.
  class GenericWriteStream {
    public:
      typedef char Ch;    //!< Character type. Only support char.

      //! Constructor
      /*!
        \param os Output stream.
        */
      GenericWriteStream(std::ostream& os) : os_(os) {
      }

      void Put(char c) {
        os_.put(c);
      }

      void PutN(char c, size_t n) {
        for (size_t i = 0; i < n; ++i) {
          Put(c);
        }
      }

      void Flush() {
        os_.flush();
      }

      size_t Tell() const {
        return (int)os_.tellp();
      }

      // Not implemented
      char Peek() const    { RAPIDJSON_ASSERT(false); }
      char Take()          { RAPIDJSON_ASSERT(false); }
      char* PutBegin()     { RAPIDJSON_ASSERT(false); return 0; }
      size_t PutEnd(char*) { RAPIDJSON_ASSERT(false); return 0; }

    private:
      std::ostream& os_;
  };

  template<>
    inline void PutN(GenericWriteStream& stream, char c, size_t n) {
      stream.PutN(c, n);
    }

} // namespace rapidjson

#endif // RAPIDJSON_GENERICSTREAM_H_
