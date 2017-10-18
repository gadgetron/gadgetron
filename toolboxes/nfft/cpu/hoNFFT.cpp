/** hoNFFT.cpp */

#include "hoNFFT.h"

#include "hoNDFFT.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_utils.h"

#include "vector_td_utilities.h"
#include "vector_td_operators.h"
#include "vector_td_io.h"

#include <math.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <boost/make_shared.hpp>

using namespace std;

namespace Gadgetron{
	
	template<class Real, unsigned int D>
	hoNFFT_plan<Real, D>::hoNFFT_plan(){
		throw std::runtime_error("Default constructor is not available");
	}

	template<class Real, unsigned int D>
	hoNFFT_plan<Real, D>::hoNFFT_plan(
		typename uint64d<D>::Type n,
		Real osf,
		Real wg
	){
		if(osf < Real(1.0) || osf > Real(2.0))
			throw std::runtime_error("Oversampling factor must be between 1 and 2");

		if(wg < Real(1.0) || wg > Real(10.0))
			throw std::runtime_error("Kernel width must be between 1 and 10");

		for(size_t i = 0; i < D; i++)
			if(n[i] < 0) 
				throw std::runtime_error("Matrix size must be positive");
		
		auto v = n[0];
		for(size_t i = 0; i < D; i++)
			if(n[i] != v)
				throw std::runtime_error("Matrix dimensions must be equal");
	
		this->n = n;
		this->osf = osf;
		this->wg = wg;
	}

	template<class Real, unsigned int D>
	hoNFFT_plan<Real, D>::~hoNFFT_plan(){
		// Empty destructor
	}

	template<class Real, unsigned int D>
	void hoNFFT_plan<Real, D>::preprocess(
		hoNDArray<typename reald<Real, D>::Type> k
	){
		if(k.get_number_of_elements() == 0)
			throw std::runtime_error("Empty Trajectory");

		for(auto it: k)
			for(size_t i = 0; i < D; i++)
				if(it[i] > Real(0.5) || it[i] < Real(-0.5))
				 throw std::runtime_error("Trajectory must be between [-0.5,0.5]");
		
		this->k = k;
		initialize();		
	}

	template<class Real, unsigned int D>
	void hoNFFT_plan<Real, D>::compute(
		hoNDArray<complext<Real>> &d,
		hoNDArray<complext<Real>> &m,
		hoNDArray<Real> w,
		NFFT_comp_mode mode
	){
		if(d.get_number_of_elements() == 0)
			throw std::runtime_error("Empty data");

		if(m.get_number_of_elements() == 0)
			throw std::runtime_error("Empty gridding matrix");

		switch(mode){
			case NFFT_FORWARDS_C2NC:{
				deapodize(d, true);
				fft(d, NFFT_FORWARDS);
				convolve(d, m, NFFT_CONV_C2NC);
				
				if(w.get_number_of_elements() != 0){
					if(m.get_number_of_elements() != w.get_number_of_elements())
						throw std::runtime_error("Incompatible dimensions");

					m /= w;
				}
				break;
			}
			case NFFT_FORWARDS_NC2C:{
				if(w.get_number_of_elements() != 0){
					if(w.get_number_of_elements() != d.get_number_of_elements())
						throw std::runtime_error("Incompitalbe dimensions");

					d *= w;
				}
				
				convolve(d, m, NFFT_CONV_NC2C);
				fft(m, NFFT_FORWARDS);
				deapodize(m, true);

				break;
			}
			case NFFT_BACKWARDS_NC2C:{
				if(w.get_number_of_elements() != 0){
					if(w.get_number_of_elements() != d.get_number_of_elements())
						throw std::runtime_error("Incompatible dimensions");
					
					d *= w;
				}
					
				convolve(d, m, NFFT_CONV_NC2C);
				fft(m, NFFT_BACKWARDS);
				deapodize(m);

				break;
			}
			case NFFT_BACKWARDS_C2NC:{
				deapodize(d, true);
				fft(d, NFFT_BACKWARDS);
				convolve(d, m, NFFT_CONV_C2NC);
				
				if(w.get_number_of_elements() != 0){
					if(w.get_number_of_elements() != m.get_number_of_elements())
						throw std::runtime_error("Incompatible dimensions");

					m *= w;
				}
				break;
			}
		};
	}

	template<class Real, unsigned int D>
	void hoNFFT_plan<Real, D>::mult_MH_M(
		hoNDArray<complext<Real>> &in,
		hoNDArray<complext<Real>> &out
	){
		hoNDArray<Real> w((size_t)0);
		hoNDArray<complext<Real>> tmp;
		tmp.create(n[0]*osf,n[1]*osf);
		compute(in, tmp, w, NFFT_BACKWARDS_NC2C);
		compute(tmp, out, w, NFFT_FORWARDS_C2NC);
	}

	template<class Real, unsigned int D>
	void hoNFFT_plan<Real, D>::convolve(
		hoNDArray<complext<Real>> &d,
		hoNDArray<complext<Real>> &m,
		NFFT_conv_mode mode
	){
		if(mode == NFFT_CONV_NC2C)
			convolve_NFFT_NC2C(d, m);
		else
			convolve_NFFT_C2NC(d, m);
	}

	template<class Real, unsigned int D>
	void hoNFFT_plan<Real, D>::fft(
		hoNDArray<complext<Real>> &d,
		NFFT_fft_mode mode
	){
		if(mode == NFFT_FORWARDS)
			hoNDFFT<Real>::instance()->fft(&d);
		else
			hoNDFFT<Real>::instance()->ifft(&d);
	}

	template<class Real, unsigned int D>
	void hoNFFT_plan<Real, D>::deapodize(
		hoNDArray<complext<Real>> &d,
		bool fourierDomain
	){
		if(fourierDomain){
			if(da.get_number_of_elements() != d.get_number_of_elements())
				throw std::runtime_error("Incompatiblef deapodization dimensions");
			
			d *= da;
		}else{
			if(da.get_number_of_elements() != d.get_number_of_elements())
				throw std::runtime_error("Incompatible deapodization dimensions");

			d /= da;
		}
	}

	template<class Real, unsigned int D>
	void hoNFFT_plan<Real, D>::initialize(){
		kw = wg/osf;
		kosf = std::floor(0.91/(osf*1e-3));
		kwidth = osf*kw/2;

		Real tmp = kw*(osf-0.5);
		beta = M_PI*std::sqrt(tmp*tmp-0.8);

		p.create(kosf*kwidth+1);
		for(size_t i = 0; i < kosf*kwidth+1; i++){
			Real om = Real(i)/Real(kosf*kwidth);
			p[i] = bessi0(beta*std::sqrt(1-om*om));
		}
		Real pConst = p[0];
		for(auto it = p.begin(); it != p.end(); it++)
			*it /= pConst;
		p[kosf*kwidth] = 0;
		
		// Need to fix to allow for flexibility in dimensions
		hoNDArray<Real> dax(osf*n[0]);
		for(int i = 0; i < osf*n[0]; i++){
			Real x = (i-osf*n[0]/2)/n[0];
			Real tmp = M_PI*M_PI*kw*kw*x*x-beta*beta;
			auto sqa = std::sqrt(complex<Real>(tmp, 0));
			dax[i] = (std::sin(sqa)/sqa).real();
		}
		auto daxConst = dax[osf*n[0]/2-1];
		for(auto it = dax.begin(); it != dax.end(); it++)
			*it /= daxConst;

		switch(D){
			case 1:{
				da.create(osf*n[0]);
				std::copy(dax.begin(), dax.end(), da.begin());
				nx.create(k.get_number_of_elements());
				for(size_t i = 0; i < k.get_number_of_elements(); i++)
					nx[i] = (n[0]*osf/2)+osf*n[0]*k[i][0];
				break;
			}
			case 2:{
				da.create(osf*n[0], osf*n[1]);
				for(size_t i = 0; i < osf*n[0]; i++)
					for(size_t j = 0; j < osf*n[1]; j++)
						da[i+j*n[1]*osf] = dax[i]*dax[j];
				nx.create(k.get_number_of_elements());
				ny.create(k.get_number_of_elements());
				for(size_t i = 0; i < k.get_number_of_elements(); i++){
					nx[i] = (n[0]*osf/2)+osf*n[0]*k[i][0];
					ny[i] = (n[1]*osf/2)+osf*n[1]*k[i][1];
				}
				break;
			}
			case 3:{
				da.create(osf*n[0], osf*n[1], osf*n[2]);
				for(size_t i = 0; i < osf*n[0]; i++)
					for(size_t j = 0; j < osf*n[1]; j++)
						for(size_t k = 0; k < osf*n[3]; k++)
							da[i+n[1]*j*osf+n[2]*k*osf] = dax[i]*dax[j]*dax[k];
				nx.create(k.get_number_of_elements());
				ny.create(k.get_number_of_elements());
				nz.create(k.get_number_of_elements());
				for(size_t i = 0; i < k.get_number_of_elements(); i++){
					nx[i] = (n[0]*osf/2)+osf*n[0]*k[i][0];
					ny[i] = (n[1]*osf/2)+osf*n[1]*k[i][1];
					nz[i] = (n[2]*osf/2)+osf*n[2]*k[i][2];
				}
				break;
				
			}
		}
	}

	template<class Real, unsigned int D>
	void hoNFFT_plan<Real, D>::convolve_NFFT_C2NC(
		hoNDArray<complext<Real>> &m,
		hoNDArray<complext<Real>> &d
	){
		switch(D){
			case 1:{
				m.fill(0);
				for(size_t i = 0; i < k.get_number_of_elements(); i++){
					for(int lx = -kwidth; lx < kwidth+1; lx++){
						Real nxt = std::round(nx[i]+lx);
						Real kkx = std::min(
							std::round(kosf*std::abs(nx[i]-nxt)),
							std::floor(kosf*kwidth)
						);
						Real kwx = p[kkx];
						nxt = std::max(nxt, Real(0)); nxt = std::min(nxt, osf*n[0]-1);
						d[i] += m[(size_t)nxt]*kwx;
					}
				}
				break;				
			}
			case 2:{
				d.fill(0);
				for(size_t i = 0; i < k.get_number_of_elements(); i++){
					for(int lx = -kwidth; lx < kwidth+1; lx++){
						for(int ly = -kwidth; ly < kwidth+1; ly++){
							Real nxt = std::round(nx[i]+lx);
							Real nyt = std::round(ny[i]+ly);

							Real kkx = std::min(
								std::round(kosf*std::abs(nx[i]-nxt)),
								std::floor(kosf*kwidth)
							);
							Real kky = std::min(
								std::round(kosf*std::abs(ny[i]-nyt)),
								std::floor(kosf*kwidth)
							);
							Real kwx = p[kkx]; Real kwy = p[kky];

							nxt = std::max(nxt, Real(0)); nxt = std::min(nxt, osf*n[0]-1);
							nyt = std::max(nyt, Real(0)); nyt = std::min(nyt, osf*n[1]-1);

							d[i] += m[(size_t)(nxt+nyt*osf*n[1])]*kwx*kwy;
						}
					}
				}
				break;
			}
			case 3:{
				m.fill(0);
				for(size_t i = 0; i < k.get_number_of_elements(); i++){
					for(int lx = -kwidth; lx < kwidth+1; lx++){
						for(int ly = -kwidth; ly < kwidth+1; ly++){
							for(int lz = -kwidth; lz = kwidth+1; lz++){
								Real nxt = std::round(nx[i]+lx);
								Real nyt = std::round(ny[i]+ly);
								Real nzt = std::round(nz[i]+lz);

								Real kkx = std::min(
									std::round(kosf*std::abs(nx[i]-nxt)),
									std::floor(kosf*kwidth)
								);
								Real kky = std::min(
									std::round(kosf*std::abs(ny[i]-nyt)),
									std::floor(kosf*kwidth)
								);
								Real kkz = std::min(
									std::round(kosf*std::abs(nz[i]-nzt)),
									std::floor(kosf*kwidth)
								);
								Real kwx = p[kkx];
								Real kwy = p[kky];
								Real kwz = p[kkz];

								nxt = std::max(nxt, Real(0));
								nxt = std::min(nxt, osf*n[0]-1);

								nyt = std::max(nxt, Real(0));
								nyt = std::min(nyt, osf*n[1]-1);

								nzt = std::max(nzt, Real(0));
								nzt = std::min(nzt, osf*n[2]-1);

								d[i] += m[(size_t)(nxt+nyt*osf*n[1]+nzt*osf*n[2])]*kwx*kwy*kwz;
							}
						}
					}
				}
				break;	
			}
		}
	}

	template<class Real, unsigned int D>
	void hoNFFT_plan<Real, D>::convolve_NFFT_NC2C(
		hoNDArray<complext<Real>> &d,
		hoNDArray<complext<Real>> &m
	){
		switch(D){
			case 1:{
				m.fill(0);
				for(size_t i = 0; i < k.get_number_of_elements(); i++){
					complext<Real> dw = d[i];
					for(int lx = -kwidth; lx < kwidth+1; lx++){
						Real nxt = std::round(nx[i]+lx);
						Real kkx = std::min(
							std::round(kosf*std::abs(nx[i]-nxt)),
							std::floor(kosf*kwidth)
						);
						Real kwx = p[kkx];
						nxt = std::max(nxt, Real(0)); nxt = std::min(nxt, osf*n[0]-1);
						m[(size_t)nxt] += dw*kwx;
					}
				}

				m[0] = 0;
				m[m.get_number_of_elements()] = 0;
				break;
			}
			case 2:{
				m.fill(0);
				for(size_t i = 0; i < k.get_number_of_elements(); i++){
					complext<Real> dw = d[i];
					for(int lx = -kwidth; lx < kwidth+1; lx++){
						for(int ly = -kwidth; ly < kwidth+1; ly++){
							Real nxt = std::round(nx[i]+lx);
							Real nyt = std::round(ny[i]+ly);

							Real kkx = std::min(
								std::round(kosf*std::abs(nx[i]-nxt)),
								std::floor(kosf*kwidth)
							);
							Real kky = std::min(
								std::round(kosf*std::abs(ny[i]-nyt)),
								std::floor(kosf*kwidth)
							);
							Real kwx = p[kkx]; Real kwy = p[kky];

							nxt = std::max(nxt, Real(0)); nxt = std::min(nxt, osf*n[0]-1);
							nyt = std::max(nyt, Real(0)); nyt = std::min(nyt, osf*n[1]-1);

							m[(size_t)(nxt+nyt*osf*n[1])] += dw*kwx*kwy;
						}
					}

					for(size_t i = 0; i < n[0]*osf; i++){
						m[i] = 0;
						m[n[0]*osf+i] = 0;
						m[n[0]*osf*(n[0]*osf-1)+i] = 0;
						m[n[0]*osf*i+(n[0]*osf-1)] = 0;
					}
				}
				break;
			}
			case 3:{
				m.fill(0);
				for(size_t i = 0; i < k.get_number_of_elements(); i++){
					complext<Real> dw = d[i];
					for(int lx = -kwidth; lx < kwidth+1; lx++){
						for(int ly = -kwidth; ly < kwidth+1; ly++){
							for(int lz = -kwidth; lz = kwidth+1; lz++){
								Real nxt = std::round(nx[i]+lx);
								Real nyt = std::round(ny[i]+ly);
								Real nzt = std::round(nz[i]+lz);

								Real kkx = std::min(
									std::round(kosf*std::abs(nx[i]-nxt)),
									std::floor(kosf*kwidth)
								);
								Real kky = std::min(
									std::round(kosf*std::abs(ny[i]-nyt)),
									std::floor(kosf*kwidth)
								);
								Real kkz = std::min(
									std::round(kosf*std::abs(nz[i]-nzt)),
									std::floor(kosf*kwidth)
								);
								Real kwx = p[kkx];
								Real kwy = p[kky];
								Real kwz = p[kkz];

								nxt = std::max(nxt, Real(0));
								nxt = std::min(nxt, osf*n[0]-1);

								nyt = std::max(nxt, Real(0));
								nyt = std::min(nyt, osf*n[1]-1);

								nzt = std::max(nzt, Real(0));
								nzt = std::min(nzt, osf*n[2]-1);

								m[(size_t)(nxt+nyt*osf*n[1]+nzt*osf*n[2])] +=
									dw*kwx*kwy*kwz;
							}
						}
					}
				}

				for(size_t i = 0; i < n[0]*osf; i++){
					m[i] = 0;
					m[n[0]*osf+i] = 0;
					m[n[0]*osf*(n[0]*osf-1)+i] = 0;
					m[n[0]*osf*i+(n[0]*osf-1)] = 0;
					// Need to add two more
				}
				break;
			}
		}
	}

	template<class Real, unsigned int D>
	Real hoNFFT_plan<Real, D>::bessi0(Real x){
		Real denominator;
		Real numerator;
		Real z;
		if (x == 0.0) {
  		return 1.0;
  	} else {
  		z = x * x;
  	  numerator = (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* 
  	  	(z* 0.210580722890567e-22  + 0.380715242345326e-19 ) +
  	   	0.479440257548300e-16) + 0.435125971262668e-13 ) +
  	    0.300931127112960e-10) + 0.160224679395361e-7  ) +
  	    0.654858370096785e-5)  + 0.202591084143397e-2  ) +
  	    0.463076284721000e0)   + 0.754337328948189e2   ) +
  	    0.830792541809429e4)   + 0.571661130563785e6   ) +
  	    0.216415572361227e8)   + 0.356644482244025e9   ) +
  	   	0.144048298227235e10);

  	  denominator = (z*(z*(z-0.307646912682801e4)+
  	    0.347626332405882e7)-0.144048298227235e10);
		}

  	return -numerator/denominator;
	}
}

template class EXPORTCPUNFFT Gadgetron::hoNFFT_plan<float, 1>;
template class EXPORTCPUNFFT Gadgetron::hoNFFT_plan<float, 2>;
template class EXPORTCPUNFFT Gadgetron::hoNFFT_plan<float, 3>;

template class EXPORTCPUNFFT Gadgetron::hoNFFT_plan<double, 1>;
template class EXPORTCPUNFFT Gadgetron::hoNFFT_plan<double, 2>;
template class EXPORTCPUNFFT Gadgetron::hoNFFT_plan<double, 3>;

