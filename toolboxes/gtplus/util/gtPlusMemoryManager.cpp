/** \file       gtPlusMemoryManager.cpp
    \brief      Implement a simple memory manager for GtPlus toolbox
    \author     Hui Xue
*/

#include <gtPlusMemoryManager.h>
#include <cstring>
#include "log.h"

namespace Gadgetron { namespace gtPlus {

gtPlusMemoryManager::gtPlusMemoryManager(size_t aligned_bytes, size_t preallocated_bytes) : aligned_bytes_(aligned_bytes)
{
    try
    {
        memory_.reserve(1024);
        memObjList_.reserve(1024);
        GADGET_CHECK_THROW(increase(preallocated_bytes));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in gtPlusMemoryManager::gtPlusMemoryManager(aligned_bytes, preallocated_bytes) : " << preallocated_bytes/1024/1024 << " MegaBytes ... ");
    }
}

gtPlusMemoryManager::~gtPlusMemoryManager()
{
    try
    {
        // release all chunks hold by manager
        size_t num = memory_.size();
        for ( size_t ii=0; ii<num; ii++ )
        {
            _deallocate_memory(memory_[ii].mem_chunk_ptr_);
            memory_[ii].len_chunk_bytes_ = 0;
        }

        // clear the list
        allocated_list_.clear();
        free_list_.clear();
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in gtPlusMemoryManager::~gtPlusMemoryManager() ... ");
    }
}

void* gtPlusMemoryManager::allocate(size_t size)
{
    void* ptr = NULL;
    mutex_.lock();
    ptr = allocateImpl(size);
    mutex_.unlock();
    return ptr;
}

void* gtPlusMemoryManager::allocateImpl(size_t size)
{
    try
    {
        if ( size == 0 )
        {
            return NULL;
        }

        // go through the free list and find a block big enough
        MemoryListType::iterator iter = free_list_.begin();
        for ( ; iter!=free_list_.end(); iter++ )
        {
            if ( iter->second.len_bytes_ >= size )
            {
                break;
            }
        }

        if ( iter == free_list_.end() )
        {
            this->defragmentImpl();

            iter = free_list_.begin();
            for ( ; iter!=free_list_.end(); iter++ )
            {
                if ( iter->second.len_bytes_ >= size )
                {
                    break;
                }
            }

            if ( iter == free_list_.end() )
            {
                // increase the managed buffer
                if ( !increase(size) ) return NULL;
                // if ( !increase(size) ) return NULL;

                // allocate the last chunk
                return this->allocateChunkAsUsed(memory_.size()-1);
            }
        }

        if ( iter != free_list_.end() )
        {
            gtPlusMemoryObj obj = iter->second;

            gtPlusMemoryObj allocateObj = obj;
            allocateObj.len_bytes_ = size;
            if ( allocateObj.len_bytes_%aligned_bytes_ != 0 )
            {
                allocateObj.len_bytes_ = aligned_bytes_ - (allocateObj.len_bytes_%aligned_bytes_);
            }

            gtPlusMemoryObj freeObj = obj;
            freeObj.len_bytes_ = obj.len_bytes_-allocateObj.len_bytes_;
            freeObj.mem_ptr_ = (void*)((char*)obj.mem_ptr_+allocateObj.len_bytes_);

            // modify the current free obj
            free_list_.erase(iter);
            free_list_[freeObj.mem_ptr_] = freeObj;

            // insert the allocated mem into allocated list
            allocated_list_[allocateObj.mem_ptr_] = allocateObj;

            return allocateObj.mem_ptr_;
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Error happened in gtPlusMemoryManager::allocate(size_t size) : " << size);
        return NULL;
    }

    return NULL;
}

void gtPlusMemoryManager::free(void* raw_memory)
{
    mutex_.lock();
    freeImpl(raw_memory);
    mutex_.unlock();
}

void gtPlusMemoryManager::freeImpl(void* raw_memory)
{
    try
    {
        MemoryListType::iterator iter;
        iter = allocated_list_.find(raw_memory);
        if(iter == allocated_list_.end()) return;

        free_list_[iter->first] = iter->second;
        allocated_list_.erase(iter);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Error happened in gtPlusMemoryManager::freeImpl(raw_memory) ... ");
    }
}

void gtPlusMemoryManager::defragment()
{
    mutex_.lock();
    defragmentImpl();
    mutex_.unlock();
}

void gtPlusMemoryManager::defragmentImpl()
{
    try
    {
        GADGET_MSG_DEPRECATED("-----> gtPlusMemoryManager::defragmentImpl() ... ");
        size_t N = free_list_.size();
        memObjList_.resize(N);

        size_t ii=0;
        MemoryListType::iterator iter = free_list_.begin();
        for ( ; iter!=free_list_.end(); iter++ )
        {
            memObjList_[ii++] = iter->second;
        }

        std::sort(memObjList_.begin(), memObjList_.end(), MemoryObjCompare() );

        size_t numChunk = memory_.size();

        size_t n, jj;

        for ( n=0; n<numChunk; n++ )
        {
            size_t start(0), end(0);
            N = memObjList_.size();
            bool rangeFound = false;

            for ( ii=0; ii<N; ii++ )
            {
                if( memObjList_[ii].chunk_id_ == n )
                {
                    start = ii;
                    for ( jj=start+1; jj<N; jj++ )
                    {
                        if( memObjList_[jj].chunk_id_ > n )
                        {
                            end = jj-1;
                            rangeFound = true;
                            break;
                        }
                    }

                    if( rangeFound )
                    {
                        break;
                    }
                }
            }

            if ( end > start )
            {
                while ( true )
                {
                    for ( ii=start; ii<end; ii++ )
                    {
                        if ( ((char*)memObjList_[ii].mem_ptr_+memObjList_[ii].len_bytes_) >= memObjList_[ii+1].mem_ptr_ )
                        {
                            // combine ii and ii+1 
                            gtPlusMemoryObj obj(memObjList_[ii]);
                            obj.len_bytes_ = (char*)(memObjList_[ii+1].mem_ptr_)+memObjList_[ii+1].len_bytes_ - (char*)memObjList_[ii].mem_ptr_;

                            memObjList_[ii] = obj;
                            memObjList_.erase(memObjList_.begin()+ii+1);

                            end--;

                            ii=start;
                            break;
                        }
                    }

                    if ( ii == end )
                    {
                        break;
                    }
                }
            }
        }

        free_list_.clear();

        N = memObjList_.size();
        for ( ii=0; ii<N; ii++ )
        {
            free_list_[memObjList_[ii].mem_ptr_] = memObjList_[ii];
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Error happened in gtPlusMemoryManager::defragmentImpl() ... ");
    }
}

void* gtPlusMemoryManager::allocateChunkAsUsed(size_t chunkId)
{
    gtPlusMemoryChunkObj mem = memory_[chunkId];

    // add the chunk to the allocated list
    gtPlusMemoryObj obj;
    obj.chunk_id_ = chunkId;
    obj.mem_ptr_ = mem.mem_chunk_ptr_;
    obj.len_bytes_ = mem.len_chunk_bytes_;

    allocated_list_[obj.mem_ptr_] = obj;

    // remove the chunk from the free list
    MemoryListType::iterator it;
    it = free_list_.find(obj.mem_ptr_);
    if(it == free_list_.end()) return NULL;
    free_list_.erase(it);

    return mem.mem_chunk_ptr_;
}

bool gtPlusMemoryManager::increase(size_t added_bytes)
{
    try
    {
        GADGET_MSG_DEPRECATED("-----> gtPlusMemoryManager::increase() : " << added_bytes/1024/1024 << " MegaBytes ");

        void* ptr;
        _allocate_memory(added_bytes, ptr);
        if ( ptr==NULL ) return false;

        // create on chuck
        gtPlusMemoryChunkObj chunk;
        chunk.len_chunk_bytes_ = added_bytes;
        chunk.mem_chunk_ptr_ = ptr;

        // insert into the lists
        memory_.push_back(chunk);

        // insert into the free list
        gtPlusMemoryObj mem;
        mem.chunk_id_ = memory_.size()-1;
        mem.mem_ptr_ = ptr;
        mem.len_bytes_ = added_bytes;

        free_list_[ptr] = mem;
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in gtPlusMemoryManager::increase(size_t added_bytes) : " << added_bytes/1024/1024 << " MegaBytes ...");
        return false;
    }

    return true;
}

void gtPlusMemoryManager::_allocate_memory( size_t size, void*& data )
{
    data = calloc(size, 1);
}

void gtPlusMemoryManager::_deallocate_memory( void* data )
{
    free(data);
}

size_t gtPlusMemoryManager::totalFreeMemory() const
{
    size_t memSize = 0;
    MemoryListType::const_iterator iter = free_list_.begin();
    for ( ; iter!=free_list_.end(); iter++ )
    {
        memSize += iter->second.len_bytes_;
    }
    return memSize;
}

size_t gtPlusMemoryManager::maxFreeMemoryChunkSize() const
{
    size_t maxChunkSize = 0;
    MemoryListType::const_iterator iter = free_list_.begin();
    for ( ; iter!=free_list_.end(); iter++ )
    {
        if ( iter->second.len_bytes_ > maxChunkSize ) maxChunkSize = iter->second.len_bytes_;
    }
    return maxChunkSize;
}

void gtPlusMemoryManager::printInfo(std::ostream& os)
{
    using namespace std;
    os << "-------------- GTPlus ISMRMRD Recon Memory Manager -------------" << endl;
    os << "Implementation of a simple memory manager for large chunk memory management" << endl;
    os << "Managed chunk : " << memory_.size() << endl;
    size_t ii;
    for ( ii=0; ii<memory_.size(); ii++ )
    {
        os << "--> Chunk " << ii << " - " << memory_[ii].len_chunk_bytes_/1024 << " kiloBytes <--" << endl;
    }
    os << "----------------------------------" << endl;
    os << "Allocated memory pieces  : " << allocated_list_.size() << endl;
    MemoryListType::iterator iter = allocated_list_.begin();
    ii=0;
    for ( ; iter!=allocated_list_.end(); iter++ )
    {
        os << "--> Allocated " << ii++ << " - " << iter->second.len_bytes_/1024 << " kiloBytes <--" << endl;
    }
    os << "----------------------------------" << endl;
    os << "Free memory pieces  : " << free_list_.size() << endl;
    iter = free_list_.begin();
    ii=0;
    for ( ; iter!=free_list_.end(); iter++ )
    {
        os << "--> Free " << ii++ << " - " << iter->second.len_bytes_/1024 << " kiloBytes <--" << endl;
    }
    os << "----------------------------------------------------------------" << endl;
}

}}
