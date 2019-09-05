#include "cuNDArray_operators.h"
#include "complext.h"

#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <complex>

namespace Gadgetron{

  // Private utility to verify array dimensions. 
  // It "replaces" NDArray::dimensions_equal() to support batch mode.
  // There is an identical function for all array instances (currently hoNDArray, cuNDArray, hoCuNDAraay)
  // !!! Remember to fix any bugs in all versions !!!
  //
  template<class T,class S> static bool compatible_dimensions( const cuNDArray<T> &x, const cuNDArray<S> &y )
  {
    return ((x.get_number_of_elements()%y.get_number_of_elements())==0);
  }

  template<typename T>
  class cuNDA_modulus : public thrust::unary_function<T,T>
  {
  public:
    cuNDA_modulus(int x):mod(x) {};
    __host__ __device__ T operator()(const T &y) const {return y%mod;}
  private:
    const int mod;
  };

  //
  // This transform support batch mode when the number of elements in x is a multiple of the number of elements in y
  //
  template<class T,class S,class F>  
  static void equals_transform(cuNDArray<T> &x, const cuNDArray<S> &y){
    if (x.dimensions_equal(&y)){
      thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), F());
    } else if (compatible_dimensions(x,y))
      {
        if (y.get_number_of_elements() < x.get_number_of_elements()) {
          typedef thrust::transform_iterator<cuNDA_modulus<int>, thrust::counting_iterator<int>, int> transform_it;
          transform_it indices = thrust::make_transform_iterator(thrust::make_counting_iterator(0),
                                                                 cuNDA_modulus<int>(y.get_number_of_elements()));
          thrust::permutation_iterator<thrust::device_ptr<S>, transform_it> p = thrust::make_permutation_iterator(
                  y.begin(), indices);
          thrust::transform(x.begin(), x.end(), p, x.begin(), F());
        } else {
          thrust::transform(x.begin(),x.end(),y.begin(),x.begin(),F());
        }

      } else {
      throw std::runtime_error("The provided cuNDArrays have incompatible dimensions for operator {+=,-=,*=,/=}");
    }
  }





  template<typename T>
  struct cuNDA_plus : public thrust::binary_function<complext<T>, T, complext<T> >
  {
    __device__ complext<T> operator()(const complext<T> &x, const T &y) const {return x+y;}
  };

  template<typename T>
  struct cuNDA_minus : public thrust::binary_function<complext<T>, T, complext<T> >
  {
    __device__ complext<T> operator()(const complext<T> &x, const T &y) const {return x-y;}
  };

  template<typename T>
  struct cuNDA_multiply : public thrust::binary_function<complext<T>, T, complext<T> >
  {
    __device__ complext<T> operator()(const complext<T> &x, const T &y) const {return x*y;}
  };

  template<typename T>
  struct cuNDA_divide : public thrust::binary_function<complext<T>, T, complext<T> >
  {
    __device__ complext<T> operator()(const complext<T> &x, const T &y) const {return x/y;}
  };

  template<class T> cuNDArray<typename std::enable_if<enable_operator<T>(), T >::type > & operator+= (cuNDArray<T> &x, const  cuNDArray<T> &y){
    equals_transform< T,T,thrust::plus<T> >(x,y);
    return x;
  }

  template<class T> cuNDArray<typename std::enable_if<enable_operator<T>(), T >::type > & operator+= (cuNDArray<T> &x , T y){
    thrust::constant_iterator<T> iter(y);
    thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::plus<T>());
    return x;
  }

  template<class T> cuNDArray<complext<typename std::enable_if<enable_operator<T>(), T >::type > >& operator+= (cuNDArray< complext<T> > &x , const cuNDArray<T> &y){
    equals_transform< complext<T>,T,cuNDA_plus<T> >(x,y);
    return x;
  }

  template<class T> cuNDArray<complext<typename std::enable_if<enable_operator<T>(), T >::type > >& operator+= (cuNDArray<complext<T> > &x , T y){
    thrust::constant_iterator<T> iter(y);
    thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_plus<T>());
    return x;
  }

  template<class T> cuNDArray<typename std::enable_if<enable_operator<T>(), T >::type >& operator-= (cuNDArray<T> & x , const cuNDArray<T> & y){
    equals_transform< T,T,thrust::minus<T> >(x,y);
    return x;
  }

  template<class T> cuNDArray<typename std::enable_if<enable_operator<T>(), T >::type >& operator-= (cuNDArray<T> &x , T y){
    thrust::constant_iterator<T> iter(y);
    thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::minus<T>());
    return x;
  }

  template<class T> cuNDArray<complext<typename std::enable_if<enable_operator<T>(), T >::type > >& operator-= (cuNDArray< complext<T> > &x , const cuNDArray<T> &y){
    equals_transform< complext<T>,T,cuNDA_minus<T> >(x,y);
    return x;
  }

  template<class T> cuNDArray<complext<typename std::enable_if<enable_operator<T>(), T >::type > >& operator-= (cuNDArray<complext<T> > &x , T y){
    thrust::constant_iterator<T> iter(y);
    thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_minus<T>());
    return x;
  }

  template<class T> cuNDArray<typename std::enable_if<enable_operator<T>(), T >::type >& operator*= (cuNDArray<T> &x , const cuNDArray<T> &y){
    equals_transform< T,T,thrust::multiplies<T> >(x,y);
    return x;
  }

  template<class T> cuNDArray<typename std::enable_if<enable_operator<T>(), T >::type >& operator*= (cuNDArray<T> &x , T y){
    thrust::constant_iterator<T> iter(y);
    thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::multiplies<T>());
    return x;
  }

  template<class T> cuNDArray<complext<typename std::enable_if<enable_operator<T>(), T >::type > >& operator*= (cuNDArray< complext<T> > &x , const cuNDArray<T> &y){
    equals_transform< complext<T>,T,cuNDA_multiply<T> >(x,y);
    return x;
  }

  template<class T> cuNDArray<complext<typename std::enable_if<enable_operator<T>(), T >::type > >& operator*= (cuNDArray<complext<T> > &x , T y){
    thrust::constant_iterator<T> iter(y);
    thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_multiply<T>());
    return x;
  }

  template<class T> cuNDArray<typename std::enable_if<enable_operator<T>(), T >::type >& operator/= (cuNDArray<T> &x , const cuNDArray<T> &y){
    equals_transform< T,T,thrust::divides<T> >(x,y);
    return x;
  }

  template<class T> cuNDArray<typename std::enable_if<enable_operator<T>(), T >::type >& operator/= (cuNDArray<T> &x , T y){
    thrust::constant_iterator<T> iter(y);
    thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::divides<T>());
    return x;
  }

  template<class T> cuNDArray<complext<typename std::enable_if<enable_operator<T>(), T >::type > >& operator/= (cuNDArray< complext<T> > &x , const cuNDArray<T> &y){
    equals_transform< complext<T>,T,cuNDA_divide<T> >(x,y);
    return x;
  }

  template<class T> cuNDArray<complext<typename std::enable_if<enable_operator<T>(), T >::type > >& operator/= (cuNDArray<complext<T> > &x , T y){
    thrust::constant_iterator<T> iter(y);
    thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_divide<T>());
    return x;
  }


  cuNDArray<bool>& operator&= (cuNDArray<bool> &x , cuNDArray<bool> &y){

    equals_transform< bool,bool,thrust::logical_and<bool> >(x,y);
    return x;
  }
  cuNDArray<bool>& operator|= (cuNDArray<bool> &x , cuNDArray<bool> &y){

    equals_transform< bool,bool,thrust::logical_or<bool> >(x,y);
    return x;
  }

  //
  // Instantiation
  //

  template EXPORTGPUCORE cuNDArray<float>& operator+=<float>(cuNDArray<float>&, const cuNDArray<float>&);
  template EXPORTGPUCORE cuNDArray<float>& operator+=<float>(cuNDArray<float>&, float);
  template EXPORTGPUCORE cuNDArray<float>& operator-=<float>(cuNDArray<float>&, const cuNDArray<float>&);
  template EXPORTGPUCORE cuNDArray<float>& operator-=<float>(cuNDArray<float>&, float);
  template EXPORTGPUCORE cuNDArray<float>& operator*=<float>(cuNDArray<float>&, const cuNDArray<float>&);
  template EXPORTGPUCORE cuNDArray<float>& operator*=<float>(cuNDArray<float>&, float);
  template EXPORTGPUCORE cuNDArray<float>& operator/=<float>(cuNDArray<float>&, const cuNDArray<float>&);
  template EXPORTGPUCORE cuNDArray<float>& operator/=<float>(cuNDArray<float>&, float);

  template EXPORTGPUCORE cuNDArray<double>& operator+=<double>(cuNDArray<double>&, const cuNDArray<double>&);
  template EXPORTGPUCORE cuNDArray<double>& operator+=<double>(cuNDArray<double>&, double);
  template EXPORTGPUCORE cuNDArray<double>& operator-=<double>(cuNDArray<double>&, const cuNDArray<double>&);
  template EXPORTGPUCORE cuNDArray<double>& operator-=<double>(cuNDArray<double>&, double);
  template EXPORTGPUCORE cuNDArray<double>& operator*=<double>(cuNDArray<double>&,  const cuNDArray<double>&);
  template EXPORTGPUCORE cuNDArray<double>& operator*=<double>(cuNDArray<double>&, double);
  template EXPORTGPUCORE cuNDArray<double>& operator/=<double>(cuNDArray<double>&,  const cuNDArray<double>&);
  template EXPORTGPUCORE cuNDArray<double>& operator/=<double>(cuNDArray<double>&, double);

  template EXPORTGPUCORE cuNDArray< complext<float> >& operator+=< complext<float> > 
  (cuNDArray< complext<float> >&,  const cuNDArray< complext<float> >&);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator+=< complext<float> > 
  (cuNDArray< complext<float> >&, complext<float>);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator-=< complext<float> > 
  (cuNDArray< complext<float> >&,  const cuNDArray< complext<float> >&);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator-=< complext<float> > 
  (cuNDArray< complext<float> >&, complext<float>);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator*=< complext<float> >
  (cuNDArray< complext<float> >&,  const cuNDArray< complext<float> >&);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator*=< complext<float> >
  (cuNDArray< complext<float> >&, complext<float>);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator/=< complext<float> > 
  (cuNDArray< complext<float> >&,  const cuNDArray< complext<float> >&);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator/=< complext<float> > 
  (cuNDArray< complext<float> >&, complext<float>);

  template EXPORTGPUCORE cuNDArray< complext<float> >& operator+=<float>(cuNDArray< complext<float> >&,  const cuNDArray<float>&);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator-=<float>(cuNDArray< complext<float> >&,  const cuNDArray<float>&);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator*=<float>(cuNDArray< complext<float> >&,  const cuNDArray<float>&);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator/=<float>(cuNDArray< complext<float> >&,  const cuNDArray<float>&);

  template EXPORTGPUCORE cuNDArray< complext<float> >& operator+=<float>(cuNDArray< complext<float> >&, float);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator-=<float>(cuNDArray< complext<float> >&, float);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator*=<float>(cuNDArray< complext<float> >&, float);
  template EXPORTGPUCORE cuNDArray< complext<float> >& operator/=<float>(cuNDArray< complext<float> >&, float);

  template EXPORTGPUCORE cuNDArray< complext<double> >& operator+=< complext<double> > 
  (cuNDArray< complext<double> >&,  const cuNDArray< complext<double> >&);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator+=< complext<double> > 
  (cuNDArray< complext<double> >&, complext<double>);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator-=< complext<double> > 
  (cuNDArray< complext<double> >&,  const cuNDArray< complext<double> >&);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator-=< complext<double> > 
  (cuNDArray< complext<double> >&, complext<double>);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator*=< complext<double> >
  (cuNDArray< complext<double> >&,  const cuNDArray< complext<double> >&);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator*=< complext<double> >
  (cuNDArray< complext<double> >&, complext<double>);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator/=< complext<double> > 
  (cuNDArray< complext<double> >&,  const cuNDArray< complext<double> >&);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator/=< complext<double> > 
  (cuNDArray< complext<double> >&, complext<double>);

  template EXPORTGPUCORE cuNDArray< complext<double> >& operator+=<double>(cuNDArray< complext<double> >&,  const cuNDArray<double>&);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator-=<double>(cuNDArray< complext<double> >&,  const cuNDArray<double>&);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator*=<double>(cuNDArray< complext<double> >&,  const cuNDArray<double>&);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator/=<double>(cuNDArray< complext<double> >&,  const cuNDArray<double>&);

  template EXPORTGPUCORE cuNDArray< complext<double> >& operator+=<double>(cuNDArray< complext<double> >&, double);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator-=<double>(cuNDArray< complext<double> >&, double);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator*=<double>(cuNDArray< complext<double> >&, double);
  template EXPORTGPUCORE cuNDArray< complext<double> >& operator/=<double>(cuNDArray< complext<double> >&, double);
}
