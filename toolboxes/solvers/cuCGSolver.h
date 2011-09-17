#pragma once
#include "cgSolver.h"
#include "cuNDArray.h"
#include "cuCGMatrixOperator.h"
#include "cuCGPreconditioner.h"
#include "gadgetron_export.h"

template <class REAL, class T> class EXPORTSOLVERS cuCGSolver : public cgSolver< REAL, T, cuNDArray<T> >
{
 public:

  //  cuCGSolver() : cgSolver< REAL, T, cuNDArray<T> >() { set_device(); }  
  cuCGSolver( int device=-1 ) : cgSolver< REAL, T, cuNDArray<T> >() { set_device(device); }

  virtual ~cuCGSolver() {}

  virtual bool pre_solve(cuNDArray<T>**);
  virtual bool post_solve(cuNDArray<T>**);
  virtual void solver_error( std::string err );

  virtual T solver_dot( cuNDArray<T>*, cuNDArray<T>* );
  virtual bool solver_clear( cuNDArray<T>* );
  virtual bool solver_scal( T, cuNDArray<T>* ); 
  virtual bool solver_axpy( T, cuNDArray<T>*, cuNDArray<T>* );

  virtual bool set_device(int device);

protected:
  int device_;
  int old_device_;
  cuNDArray<T> *new_rhs;
};
