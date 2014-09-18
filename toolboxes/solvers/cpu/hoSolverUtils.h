#pragma once

#include "hoNDArray.h"
#include "hoNDArray_math.h"
#include "complext.h"

#ifdef USE_OMP
#include <omp.h>
#endif

namespace Gadgetron {
template<class T> void solver_non_negativity_filter(hoNDArray<T> *xdata, hoNDArray<T> *gdata)
{
	typedef typename realType<T>::Type REAL;

	T* x = xdata->get_data_ptr();
	T* g = gdata->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
	for( int i=0; i < xdata->get_number_of_elements(); i++ )
		if( (real(x[i]) <= REAL(0)) && (real(g[i]) > 0) )
			g[i]=T(0);
}
}
