#include "cuCG.h"


template <class T> cuNDArray<T> cuCG<T>::solve(cuNDArray<T>* rhs)
{
  cuNDArray<T> ret;

  //TODO: Multiply rhs with preconditioning matrix
  //T rr_0;
  

  /*
  r = r.dot_multiply(D);

  ComplexFloat rr_0 = inner_product(r);

  ComplexFloat rr_1 = 0;
  ComplexFloat rr = 0;
  
  float rr_last = 1e10;

  ComplexFloatArray p,p_1;
  ComplexFloatArray q;

  ComplexFloatArray rho(r.get_dimensions());

  ComplexFloat alpha, beta;

  float rel_res;

  std::cout << "Iterating..." << std::endl;
  for (int it = 0; it < parm->max_iterations; it++) {
    rr_1 = rr;
    rr = inner_product(r);
    
    if (it == 0){
      p = r;
    } else {        
      beta = rr/rr_1;
      p *= beta;
      for (unsigned int i = 0; i < p.get_number_of_elements(); i++) p[i] += r[i];
    }

    p_1 = p.dot_multiply(D);
    q = mult_E_fcn(p_1, co, csm, recon_dimensions, displacementField);
    q = mult_EH_fcn(q, co, csm, recon_dimensions,
		lower_bounds,upper_bounds, image_indices, image_weights);

    ComplexFloatArray tmp_Lm = p_1.dot_multiply(theta_int);
    for (unsigned long i = 0; i < q.get_number_of_elements(); i++) {
      q[i] += std::complex<float>(lambda,0.0)*tmp_Lm[i];
    }

    q = q.dot_multiply(D);

    ComplexFloat pq = 0;
    for (unsigned int i = 0; i < q.get_number_of_elements(); i++) pq += conj(p[i])*q[i];

    alpha = rr/(pq);

    //Update current solution
    for (unsigned int i = 0; i< rho.get_number_of_elements(); i++) rho[i] += alpha*p[i];
    */

  return ret;
}

template class cuCG<float>;
template class cuCG<float2>;
