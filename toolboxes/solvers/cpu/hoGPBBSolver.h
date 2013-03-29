#pragma once

#include "hoGTBLAS.h"
#include "gpBBSolver.h"
#include "hoNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"


namespace Gadgetron{
template <class T> class hoGPBBSolver
: public gpBBSolver<hoNDArray<T> >{

public:
	typedef typename realType<T>::type REAL;
	hoGPBBSolver() : gpBBSolver<hoNDArray<T> >() { };
	virtual	~hoGPBBSolver(){};


	virtual void solver_non_negativity_filter(hoNDArray<T> *xdata,hoNDArray<T> *gdata){
		T* x = xdata->get_data_ptr();
		T* g = gdata->get_data_ptr();
		for (int i = 0; i != xdata->get_number_of_elements(); ++i)
			if ( real(x[i]) <= REAL(0) && real(g[i]) > 0) g[i]=T(0);
	}

	virtual void solver_reciprocal_clamp( hoNDArray<T>* xdata,REAL threshold){
		T* x = xdata->get_data_ptr();
		for (int i = 0; i != xdata->get_number_of_elements(); ++i){
			if (real(x[i]) < threshold) x[i]= T(0);
			else x[i]= T(1)/x[i];
		}
	}




};
}

