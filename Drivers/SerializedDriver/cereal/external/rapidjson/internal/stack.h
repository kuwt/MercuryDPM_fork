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

#ifndef RAPIDJSON_INTERNAL_STACK_H_
#define RAPIDJSON_INTERNAL_STACK_H_

namespace rapidjson {
namespace internal {

///////////////////////////////////////////////////////////////////////////////
// Stack

//! A type-unsafe stack for storing different types of data.
/*! \tparam Allocator Allocator for allocating stack memory.
*/
template <typename Allocator>
class Stack {
public:
	Stack(Allocator* allocator, size_t stack_capacity) : allocator_(allocator), own_allocator_(0), stack_(0), stack_top_(0), stack_end_(0), stack_capacity_(stack_capacity) {
		RAPIDJSON_ASSERT(stack_capacity_ > 0);
		if (!allocator_)
			own_allocator_ = allocator_ = new Allocator();
		stack_top_ = stack_ = (char*)allocator_->Malloc(stack_capacity_);
		stack_end_ = stack_ + stack_capacity_;
	}

	~Stack() {
		Allocator::Free(stack_);
		delete own_allocator_; // Only delete if it is owned by the stack
	}

	void Clear() { /*stack_top_ = 0;*/ stack_top_ = stack_; }

	template<typename T>
	T* Push(size_t count = 1) {
		 // Expand the stack if needed
		if (stack_top_ + sizeof(T) * count >= stack_end_) {
			size_t new_capacity = stack_capacity_ * 2;
			size_t size = GetSize();
			size_t new_size = GetSize() + sizeof(T) * count;
			if (new_capacity < new_size)
				new_capacity = new_size;
			stack_ = (char*)allocator_->Realloc(stack_, stack_capacity_, new_capacity);
			stack_capacity_ = new_capacity;
			stack_top_ = stack_ + size;
			stack_end_ = stack_ + stack_capacity_;
		}
		T* ret = (T*)stack_top_;
		stack_top_ += sizeof(T) * count;
		return ret;
	}

	template<typename T>
	T* Pop(size_t count) {
		RAPIDJSON_ASSERT(GetSize() >= count * sizeof(T));
		stack_top_ -= count * sizeof(T);
		return (T*)stack_top_;
	}

	template<typename T>
	T* Top() { 
		RAPIDJSON_ASSERT(GetSize() >= sizeof(T));
		return (T*)(stack_top_ - sizeof(T));
	}

	template<typename T>
	T* Bottom() { return (T*)stack_; }

	Allocator& GetAllocator() { return *allocator_; }
	size_t GetSize() const { return stack_top_ - stack_; }
	size_t GetCapacity() const { return stack_capacity_; }

private:
	Allocator* allocator_;
	Allocator* own_allocator_;
	char *stack_;
	char *stack_top_;
	char *stack_end_;
	size_t stack_capacity_;
};

} // namespace internal
} // namespace rapidjson

#endif // RAPIDJSON_STACK_H_
