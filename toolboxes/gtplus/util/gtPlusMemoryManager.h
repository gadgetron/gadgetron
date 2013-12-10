/** \file       gtPlusMemoryManager.h
    \brief      Implement a simple memory manager for GtPlus toolbox
    \author     Hui Xue
*/

#pragma once

#include "GadgetronCommon.h"
#include <iostream>
#include <typeinfo>

#include <boost/thread/mutex.hpp>
#include "GtPlusExport.h"
#include "GadgetronTimer.h"
#include "GadgetronException.h"


#include <map>
#include <vector>
#include <typeinfo>
#include <string>
#include <limits>

#include <mkl.h>

// the memory manager for large chunk allocation

namespace Gadgetron { namespace gtPlus {

struct gtPlusMemoryObj
{
    // chunk number holding this memory
    size_t chunk_id_;

    // starting address of this memory
    void* mem_ptr_;

    // memory size in bytes
    size_t len_bytes_;
};

struct gtPlusMemoryChunkObj
{
    // chunk starting address
    void* mem_chunk_ptr_;

    // memory size in bytes
    size_t len_chunk_bytes_;
};

struct MemoryObjCompare
{
    MemoryObjCompare() {}
    ~MemoryObjCompare() {}

    bool operator()(const gtPlusMemoryObj& a, const gtPlusMemoryObj& b) const
    {
        if ( a.chunk_id_ == b.chunk_id_ )
        {
            return (a.mem_ptr_ <= b.mem_ptr_);
        }
        else
        {
            return (a.chunk_id_ <= b.chunk_id_);
        }
    }
};

class EXPORTGTPLUS gtPlusMemoryManager
{
public:

    typedef std::map<void*, gtPlusMemoryObj> MemoryListType;

    gtPlusMemoryManager(size_t aligned_bytes=4, size_t preallocated_bytes=4294967296);
    virtual ~gtPlusMemoryManager();

    // allocate memory
    void* allocate(size_t size);

    // free memory
    void free(void* raw_memory);

    // increase the managed memory size
    virtual bool increase(size_t added_bytes);

    // defragment, combine the free memory obj if possible
    // whenever a memory requirement cannot be fullfilled, this funtion is called once before allocating new memory
    // user can call this function explictly to improve the efficiency
    void defragment();

    // total amount of free memory
    size_t totalFreeMemory() const;

    // maximal free memory chunk
    size_t maxFreeMemoryChunkSize() const;

    // print out the memory manager information
    void printInfo(std::ostream& os);

protected:

    // the allocated list
    MemoryListType allocated_list_;

    // the free list
    MemoryListType free_list_;

    // memory chunk hold by the manager
    std::vector<gtPlusMemoryChunkObj> memory_;

    // store the copy of memory objects
    std::vector<gtPlusMemoryObj> memObjList_;

    // aligned bytes
    size_t aligned_bytes_;

    // make sure the allocate and free are thread-safe
    boost::mutex mutex_;

    // allocate memory
    void* allocateImpl(size_t size);

    // free memory
    void freeImpl(void* raw_memory);

    // defragment
    void defragmentImpl();

    // perform memory allocation and release with system calls
    void _allocate_memory( size_t size, void*& data );
    void _deallocate_memory( void* data );

    // allocate specified chunk as allocated memory
    // return the mem_chunk_ptr_ for this chunk
    void* allocateChunkAsUsed(size_t chunkId);
};

}}
