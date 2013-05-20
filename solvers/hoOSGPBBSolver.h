#pragma once

#include "osGPBBSolver.h"
#include "hoNDArray_fileio.h"
#include <sstream>

namespace Gadgetron{

template<class ARRAY_TYPE> class hoOSGPBBSolver : public osGPBBSolver<ARRAY_TYPE>{

	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::Type REAL;
public:
	virtual ~hoOSGPBBSolver(){};
protected:
	virtual void solver_non_negativity_filter(ARRAY_TYPE *xdata,ARRAY_TYPE *gdata){
		ELEMENT_TYPE* x = xdata->get_data_ptr();
		ELEMENT_TYPE* g = gdata->get_data_ptr();
		for (int i = 0; i != xdata->get_number_of_elements(); ++i)
			if ( real(x[i]) <= REAL(0) && real(g[i]) > 0) g[i]=ELEMENT_TYPE(0);
	}

	virtual void iteration_callback(ARRAY_TYPE* x, int i){
		std::stringstream ss;
		ss << "iteration-" << i << ".real";
		write_nd_array(x,ss.str().c_str());
	}
};
/**
 * Specialization to disable class for cuNDArrays
 */
template<class T> class hoOSGPBBSolver<cuNDArray<T> > {
	virtual void hoOSGPBBSolver_only_works_for_host_arrays() =0;
};
}
