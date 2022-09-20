#pragma once
#include "complext.h"
#include "diagonalOperator.h"
#include "identityOperator.h"

#include <boost/make_shared.hpp>

namespace Gadgetron{
template <class ARRAY_TYPE> class eigenTester {


public:
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;
  eigenTester(){
	  tolerance = REAL(1e-8);
  }
  virtual ~eigenTester(){}

  ELEMENT_TYPE get_dominant_eigenvalue(){
	  boost::shared_ptr<ARRAY_TYPE> eigenVector = get_dominant_eigenvector();
	  return get_eigenvalue_from_vector(eigenVector.get());
  }


  ELEMENT_TYPE get_smallest_eigenvalue(ELEMENT_TYPE dominant_eigenvalue){
	  ELEMENT_TYPE beta = dominant_eigenvalue*2;
	  auto id_operator = boost::make_shared<identityOperator<ARRAY_TYPE>>();
	  id_operator->set_weight(-abs(beta));


	  regularization_operators_.push_back(id_operator);
	  GINFO(std::string("ID operator weight " + std::to_string(beta) + "\n").c_str());
	  ELEMENT_TYPE eig1 = get_dominant_eigenvalue();
	  regularization_operators_.pop_back();
	  return eig1+beta;


  }
  ELEMENT_TYPE get_smallest_eigenvalue(){
  	  ELEMENT_TYPE eig = get_dominant_eigenvalue();
  	  return get_smallest_eigenvalue(eig);
    }
  // Add encoding operator to solver (only one allowed)
   inline bool add_encoding_operator( boost::shared_ptr< linearOperator<ARRAY_TYPE> > op)
   {
     if( !op.get() ){
       GINFO("Error: linearSolver::add_matrix_operator : NULL operator provided\n");
       return false;
     }

     encoding_operator_ = op;

     return true;
   }
   inline void set_tolerance(REAL tolerance){
	   this->tolerance = tolerance;
   }
   // Add linear operator to solver (in addition to the encoding operator)
   inline bool add_linear_operator( boost::shared_ptr< linearOperator<ARRAY_TYPE> > op)
     {
       if( !op.get() ){
    	   GINFO("Error: linearSolver::add_matrix_operator : NULL operator provided\n");
         return false;
       }

       regularization_operators_.push_back(op);

       return true;
     }
	protected:
	 bool mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out )
	  {
	    // Basic validity checks
	    if( !in || !out ){
	      GINFO("Error: cgSolver::mult_MH_M : invalid input pointer(s)\n");
	      return false;
	    }
	    if( in->get_number_of_elements() != out->get_number_of_elements() ){
	    	GINFO("Error: cgSolver::mult_MH_M : array dimensionality mismatch\n");
	      return false;
	    }

	    // Intermediate storage
	    ARRAY_TYPE q( in->get_dimensions());
	   // Start by clearing the output
	    clear( out );

	    //Use encoding operator

	    this->encoding_operator_->mult_MH_M( in, &q, false );
	    axpy( this->encoding_operator_->get_weight(), &q, out );

	    // Iterate over regularization operators
	    for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){
		  this->regularization_operators_[i]->mult_MH_M( in, &q, false );
	      axpy( this->regularization_operators_[i]->get_weight(), &q, out );
	    }

	    return true;
	  }
	 ELEMENT_TYPE get_eigenvalue_from_vector(ARRAY_TYPE* eigenVector){
		 ARRAY_TYPE out(*eigenVector);
		 clear(&out);
		 mult_MH_M(eigenVector,&out);
		 size_t eigMax = amax(eigenVector);
		 ELEMENT_TYPE dom1 = eigenVector->at(eigMax);
		 size_t oMax = amax(&out);
		 ELEMENT_TYPE dom2 = out[oMax];
		 return dom2/dom1;

	 }

	  boost::shared_ptr<ARRAY_TYPE> get_dominant_eigenvector(){
		  GINFO_STREAM("Starting dominant eigenvector calculations " << tolerance << std::endl);
		  ELEMENT_TYPE norm = ELEMENT_TYPE(1);
		  ELEMENT_TYPE norm_old = ELEMENT_TYPE(2);

		  ARRAY_TYPE* in = new ARRAY_TYPE;
		  std::vector<size_t> image_dims = *this->encoding_operator_->get_domain_dimensions();

		  in->create(&image_dims);

		  fill(in,ELEMENT_TYPE(1));

		  ARRAY_TYPE* out = new ARRAY_TYPE;
		  out->create(&image_dims);

		  while (abs(norm-norm_old)/abs(norm)> tolerance){
			  norm_old=norm;
			  mult_MH_M(in,out);
			  GINFO_STREAM(dot(in,out) << std::endl);

			  norm = nrm2(out);

			  *out /= norm;
			  ARRAY_TYPE* tmp = in;
			  in = out;
			  out = tmp;

			  }
		  GINFO("Done\n");
		  delete in;
		  return boost::shared_ptr<ARRAY_TYPE>(out);
		}


	protected:

	  // Single encoding operator
	  boost::shared_ptr< linearOperator< ARRAY_TYPE> > encoding_operator_;
	  REAL tolerance;
	  // Vector of linear regularization operators
	  std::vector< boost::shared_ptr< linearOperator< ARRAY_TYPE> > > regularization_operators_;

};
}
