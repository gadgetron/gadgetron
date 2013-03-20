#include "hoNDArray_blas.h"
#include "complext.h"
#include <complex>

using namespace Gadgetron;

template<class T> T 
Gadgetron::dot( hoNDArray<T> *x, hoNDArray<T> *y, bool cc )
{
  arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
  arma::Col<typename stdType<T>::Type> yM = as_arma_col(y);
  typename stdType<T>::Type res = (cc) ? arma::cdot(xM,yM) : arma::dot(xM,yM);
  return *((T*)(&res));
}

template<class T> typename realType<T>::Type 
Gadgetron::asum( hoNDArray<T> *x )
{
  typedef typename realType<T>::Type realT;
  arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
  return realT(arma::norm(xM,1));
}

template<class T> typename realType<T>::Type 
Gadgetron::nrm2( hoNDArray<T> *x )
{
  typedef typename realType<T>::Type realT;
  arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
  return realT(arma::norm(xM,2));
}
  
template<class T> unsigned int 
Gadgetron::amin( hoNDArray<T> *x, T *val = 0x0 )
{
  typedef typename stdType<T>::Type aT;
  arma::Col<aT> xM = as_arma_col(x);
  unsigned int idx;
  aT min = xM.min(idx);
  if( val ) *val = *((T*)(&min));
  return idx;
}

template<class T> unsigned int 
Gadgetron::amax( hoNDArray<T> *x, T *val = 0x0 )
{
  typedef typename stdType<T>::Type aT;
  arma::Col<aT> xM = as_arma_col(x);
  unsigned int idx;
  aT max = xM.max(idx);
  if( val ) *val = *((T*)(&max));
  return idx;
}

template<class T> void 
Gadgetron::axpy( T a, hoNDArray<T> *x, hoNDArray<T> *y )
{
  typedef typename stdType<T>::Type aT;
  arma::Col<aT> xM = as_arma_col(x);
  arma::Col<aT> yM = as_arma_col(y);
  aT a2 = *((aT*)(&a));
  yM += (a2*xM);
}

//
// Instantiation
//

template EXPORTCPUCORE float dot<float>( hoNDArray<float>*, hoNDArray<float>*, bool );
template EXPORTCPUCORE float asum<float>( hoNDArray<float>* );
template EXPORTCPUCORE float nrm2<float>( hoNDArray<float>* );
template EXPORTCPUCORE unsigned int amin<float>( hoNDArray<float>*, float* );
template EXPORTCPUCORE unsigned int amax<float>( hoNDArray<float>*, float* );
template EXPORTCPUCORE void axpy<float>( float, hoNDArray<float>*, hoNDArray<float>* );

template EXPORTCPUCORE double dot<double>( hoNDArray<double>*, hoNDArray<double>*, bool );
template EXPORTCPUCORE double asum<double>( hoNDArray<double>* );
template EXPORTCPUCORE double nrm2<double>( hoNDArray<double>* );
template EXPORTCPUCORE unsigned int amin<double>( hoNDArray<double>*, double* );
template EXPORTCPUCORE unsigned int amax<double>( hoNDArray<double>*, double* );
template EXPORTCPUCORE void axpy<double>( double, hoNDArray<double>*, hoNDArray<double>* );

template EXPORTCPUCORE std::complex<float> dot< std::complex<float> >( hoNDArray< std::complex<float> >*, hoNDArray< std::complex<float> >*, bool );
template EXPORTCPUCORE float asum< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE float nrm2< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE unsigned int amin< std::complex<float> >( hoNDArray< std::complex<float> >*, std::complex<float>* );
template EXPORTCPUCORE unsigned int amax< std::complex<float> >( hoNDArray< std::complex<float> >*, std::complex<float>* );
template EXPORTCPUCORE void axpy< std::complex<float> >( std::complex<float> , hoNDArray< std::complex<float> >*, hoNDArray< std::complex<float> >* );

template EXPORTCPUCORE std::complex<double> dot< std::complex<double> >( hoNDArray< std::complex<double> >*, hoNDArray< std::complex<double> >*, bool );
template EXPORTCPUCORE double asum< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE double nrm2< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE unsigned int amin< std::complex<double> >( hoNDArray< std::complex<double> >*, std::complex<double>* );
template EXPORTCPUCORE unsigned int amax< std::complex<double> >( hoNDArray< std::complex<double> >*, std::complex<double>* );
template EXPORTCPUCORE void axpy< std::complex<double> >( std::complex<double> , hoNDArray< std::complex<double> >*, hoNDArray< std::complex<double> >* );

template EXPORTCPUCORE complext<float> dot< complext<float> >( hoNDArray< complext<float> >*, hoNDArray< complext<float> >*, bool );
template EXPORTCPUCORE float asum< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE float nrm2< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE unsigned int amin< complext<float> >( hoNDArray< complext<float> >*, complext<float>* );
template EXPORTCPUCORE unsigned int amax< complext<float> >( hoNDArray< complext<float> >*, complext<float>* );
template EXPORTCPUCORE void axpy< complext<float> >( complext<float> , hoNDArray< complext<float> >*, hoNDArray< complext<float> >* );

template EXPORTCPUCORE complext<double> dot< complext<double> >( hoNDArray< complext<double> >*, hoNDArray< complext<double> >*, bool );
template EXPORTCPUCORE double asum< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE double nrm2< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE unsigned int amin< complext<double> >( hoNDArray< complext<double> >*, complext<double>* );
template EXPORTCPUCORE unsigned int amax< complext<double> >( hoNDArray< complext<double> >*, complext<double>* );
template EXPORTCPUCORE void axpy< complext<double> >( complext<double> , hoNDArray< complext<double> >*, hoNDArray< complext<double> >* );
