#pragma once

#include <boost/smart_ptr.hpp>
#include <string>
#include <iostream>
#include "solvers_export.h"

template <class ARRAY_TYPE_IN, class ARRAY_TYPE_OUT> class solver
{
public:

  // Constructor/destructor
  solver() { output_mode_ = OUTPUT_SILENT; }
  virtual ~solver() {}
  
  // Output modes
  enum solverOutputModes { OUTPUT_SILENT = 0, OUTPUT_WARNINGS = 1, OUTPUT_VERBOSE = 2, OUTPUT_MAX = 3 };
  
  // Set/get output mode
  virtual int get_output_mode() { return output_mode_; }
  virtual void set_output_mode( int output_mode ) {
    if( !(output_mode >= OUTPUT_MAX || output_mode < 0 )) 
      output_mode_ = output_mode;
  }
  
  // Set/get starting solution/estimate for solver
  virtual void set_x0( boost::shared_ptr<ARRAY_TYPE_OUT> x0 ){ x0_ = x0; }
  virtual boost::shared_ptr<ARRAY_TYPE_OUT> get_x0(){ return x0_; }

  // Default error output
  virtual void solver_error( std::string msg ) { std::cerr << msg << std::endl; }

  // Default warning output
  virtual void solver_warning( std::string msg ) { std::cerr << msg << std::endl; }

  // Invoke solver
  virtual boost::shared_ptr<ARRAY_TYPE_OUT> solve( ARRAY_TYPE_IN* ) = 0;
 
  void* operator new(size_t bytes) { return ::new char[bytes]; }
  void* operator new(size_t s, void * p) { return p; }
  void operator delete(void *ptr) { delete[] static_cast<char*> (ptr); }

protected:
  int output_mode_;
  boost::shared_ptr<ARRAY_TYPE_OUT> x0_;
};
