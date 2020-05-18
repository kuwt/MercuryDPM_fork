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

#ifndef RAPIDJSON_STRINGBUFFER_H_
#define RAPIDJSON_STRINGBUFFER_H_

#include "rapidjson.h"
#include "internal/stack.h"

namespace rapidjson {

//! Represents an in-memory output stream.
/*!
	\tparam Encoding Encoding of the stream.
	\tparam Allocator type for allocating memory buffer.
	\implements Stream
*/
template <typename Encoding, typename Allocator = CrtAllocator>
struct GenericStringBuffer {
	typedef typename Encoding::Ch Ch;

	GenericStringBuffer(Allocator* allocator = 0, size_t capacity = kDefaultCapacity) : stack_(allocator, capacity) {}

	void Put(Ch c) { *stack_.template Push<Ch>() = c; }

	void Clear() { stack_.Clear(); }

	const char* GetString() const {
		// Push and pop a null terminator. This is safe.
		*stack_.template Push<Ch>() = '\0';
		stack_.template Pop<Ch>(1);

		return stack_.template Bottom<Ch>();
	}

	size_t Size() const { return stack_.GetSize(); }

	static const size_t kDefaultCapacity = 256;
	mutable internal::Stack<Allocator> stack_;
};

typedef GenericStringBuffer<UTF8<> > StringBuffer;

//! Implement specialized version of PutN() with memset() for better performance.
template<>
inline void PutN(GenericStringBuffer<UTF8<> >& stream, char c, size_t n) {
	memset(stream.stack_.Push<char>(n), c, n * sizeof(c));
}

} // namespace rapidjson

#endif // RAPIDJSON_STRINGBUFFER_H_
