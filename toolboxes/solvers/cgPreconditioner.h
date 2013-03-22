/** \file cgPreconditioner.h
    \brief Base class for preconditioners for the cgSolver class.
*/

#ifndef CGPRECONDITIONER_H
#define CGPRECONDITIONER_H
#pragma once

namespace Gadgetron{

  template <class ARRAY_TYPE> class cgPreconditioner
  {
  public:
    
    cgPreconditioner() {}
    virtual ~cgPreconditioner() {}
    
    virtual void apply( ARRAY_TYPE* in, ARRAY_TYPE* out) = 0;
    
    void* operator new (size_t bytes) { return ::new char[bytes]; }
    void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
    void * operator new(size_t s, void * p) { return p; }    
  };
}

#endif //CGPRECONDITIONER_H
