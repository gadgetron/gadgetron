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

template<class T> void updateF(hoNDArray<T>& data, typename realType<T>::Type alpha ,typename realType<T>::Type sigma){

  using REAL = typename realType<T>::Type;
  Gadgetron::transform(data,data,[&](auto val){
    auto tmp = val/(1+alpha*sigma);
    return tmp/std::max(REAL(1),abs(tmp));});
}




template<class T> void updateFgroup(std::vector<hoNDArray<T> >& datas, typename realType<T>::Type alpha ,typename realType<T>::Type sigma){
    using REAL = typename realType<T>::Type;

    for (int64_t i = 0; i <datas.front().size(); i++){
      REAL square = 0;
      for (int k = 0; k < datas.size(); k++){
        using namespace std;
        square += norm(datas[k][i]);
      }

      for (int k = 0; k < datas.size(); k++){
        using namespace std;
        datas[k][i] /= (1+alpha*sigma)*std::max(REAL(1),square/(1+alpha*sigma));
      }


    }
}
}

