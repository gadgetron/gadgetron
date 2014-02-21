#include "proton_kernels.h"
#include "vector_td_utilities.h"

#include "cuNDArray.h"

#include <stdio.h>

//TODO: Get rid of these defines.
#define INT_STEPS 2048
#define MAXSTEP 512
#define STEPS 3

using namespace Gadgetron;

/*template <typename T> __inline__ __host__ __device__ T sgn(T val)
{
    return copysign(T(1),val);
}


template< class T, unsigned int D > __inline__ __host__ __device__ vector_td<T,D> sgn ( const vector_td<T,D> &v1)
{
  vector_td<T,D> res;
  for(unsigned int i=0; i<D; i++ ) res.vec[i] =sgn(v1.vec[i]);
  return res;
}*/
/*
template< class T, class R, unsigned int D > __inline__ __host__ __device__ vector_td<typename vectorTDReturnType<T,R>::type,D> operator* ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
{
  vector_td<typename vectorTDReturnType<T,R>::type,D> res;
  for(unsigned int i=0; i<D; i++ )  res.vec[i]=v1.vec[i]*v2.vec[i];
  return res;
}

template< class T, class R, unsigned int D > __inline__ __host__ __device__ vector_td<typename vectorTDReturnType<T,R>::type,D> operator/ ( const vector_td<T,D> &v1, const vector_td<R,D> &v2 )
{
  vector_td<typename vectorTDReturnType<T,R>::type,D> res;
  for(unsigned int i=0; i<D; i++ )  res.vec[i]=v1.vec[i]/v2.vec[i];
  return res;
}
*/
template< class T, unsigned int D > __inline__ __host__ __device__ vector_td<T,D> remove_neg ( const vector_td<T,D> &v1)
{
  vector_td<T,D> res;
  for(unsigned int i=0; i<D; i++ ){
	  if (v1.vec[i]<0){ res.vec[i] =0;}
	  else {res.vec[i]=v1.vec[i];}
  }

  return res;
}




template <class REAL> __global__ void Gadgetron::forward_kernel(const REAL*  __restrict__ image, REAL* __restrict__ projections,
		const vector_td<REAL,3> * __restrict__ splines,  const vector_td<REAL,3> dims,
		const typename intd<3>::Type ndims, const int proj_dim, const int offset){

	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x+offset;
	if (idx < proj_dim){
		const int sid = idx*4;
		vector_td<int,3> co;


		int id,id_old;
		REAL t;
		REAL res = 0;
		//REAL length = lengths[idx];
		//REAL length = lengths[idx];
		REAL length=0;

		//Load in points to registers
		vector_td<REAL,3> a = splines[sid];
		vector_td<REAL,3> b = splines[sid+1];
		vector_td<REAL,3> c = splines[sid+2];
		vector_td<REAL,3> d = splines[sid+3];


		vector_td<REAL,3> p;
		vector_td<REAL,3> p_old=d;
		co = vector_td<int,3>((p_old+dims/2)*ndims/dims);
		co = amax(amin(co,ndims-1),0);
		id_old=co_to_idx(co,ndims);

		int steps =max(ndims)*STEPS;
		for (int i = 1; i < steps+1; i++){
			t = REAL(i)/(steps);
			p = d+t*(c+t*(b+t*a));

			//co = to_intd((p+dims/2)*ndims/dims);

			co = vector_td<int,3>((p+dims/2)*ndims/dims);
			co = amax(amin(co,ndims-1),0);
			id=co_to_idx(vector_td<int,3>(co),ndims);
			//id=co_to_idx(co,ndims);
			//REAL step_length = norm((-1.0/(steps*steps*steps)-3*t*t/steps+3*t/(steps*steps))*a+(1.0/(steps*steps)-2*t/steps)*b-c/steps);
			length += norm(p-p_old)/2;

			if(id_old != id){
				//if (min(co) >= 0 && co < ndims ) res+=image[id_old]*length;
				res+=image[id_old]*length;
				length=0;
			}

			length+= norm(p-p_old)/2;
			id_old=id;
			p_old=p;
			//co = to_intd((p+dims/2)*ndims/dims);
			//co = amax(amin(co,ndims-1),0);
		}
		projections[idx] += res;
	}

}

template <class REAL> __global__ void Gadgetron::backwards_kernel(const REAL* __restrict__ projections, REAL* __restrict__ image,
		const vector_td<REAL,3> * __restrict__ splines,  const vector_td<REAL,3> dims,
		const typename intd<3>::Type ndims, const int proj_dim, const int offset){

	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x+offset;
	if (idx < proj_dim){
		const int sid = idx*4;
		vector_td<int,3> co;

		int id,id_old;
		REAL t;

		//REAL length = lengths[idx];
		REAL length=0;

		//Load in points to registers
		vector_td<REAL,3> a = splines[sid];
		vector_td<REAL,3> b = splines[sid+1];
		vector_td<REAL,3> c = splines[sid+2];
		vector_td<REAL,3> d = splines[sid+3];
		REAL proj = projections[idx];

		vector_td<REAL,3> p;
		vector_td<REAL,3> p_old=d;
		co = vector_td<int,3>((d+dims/2)*ndims/dims);
		co = amax(amin(co,ndims-1),0);
		id_old=co_to_idx(co,ndims);

		int steps =max(ndims)*STEPS;
		for (int i = 1; i < steps; i++){
			t = REAL(i)/(steps);
			p = d+t*(c+t*(b+t*a));
			co = vector_td<int,3>((p+dims/2)*ndims/dims);
			co = amax(amin(co,ndims-1),0);
			id=co_to_idx(co,ndims);
			//REAL step_length = norm(((dt*dt*dt)+3*t*t*dt+3*t*dt*dt)*a+(dt*dt+2*t*dt)*b+c*dt);
			length += norm(p-p_old)/2;

			if(id_old != id){
				//if (min(co) >= 0 && co < ndims ) atomicAdd(&(image[id_old]),length*proj);
				atomicAdd(&(image[id_old]),length*proj);
				length=0;
			}

			length+=norm(p-p_old)/2;
			id_old=id;
			p_old=p;
			//co = to_intd((p+dims/2)*ndims/dims);
			//co = amax(amin(co,ndims-1),0);
		}

	}

}

/*
template <class REAL> __global__ void backwards_kernel(REAL* projections, REAL* image,
		vector_td<REAL,3> * splines,  const vector_td<REAL,3> dims,
		const typename uintd<3>::Type ndims, const int proj_dim){

	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < proj_dim){
		const int sid = idx*4;
			vector_td<int,3> co;
			int id;
			REAL t;
			//REAL length = lengths[idx];
			//Load in points to registers and calculate spline coefficients
			vector_td<REAL,3> a = splines[sid];
			vector_td<REAL,3> b = splines[sid+1];
			vector_td<REAL,3> c = splines[sid+2];
			vector_td<REAL,3> d = splines[sid+3];
			REAL proj = projections[idx];



			vector_td<REAL,3> p;
			vector_td<REAL,3> p_old=d;
			co = to_intd(((p_old+(dims/2))*ndims)/dims);

			vector_td<REAL,3> dir,planes,tn;

			t=0;
			//for (int j = 0; j < max(ndims)*2; j++){
			while (t < 1.0){
				dir = sgn(3*a*t*t+2*b*t+c);
				tn.vec[0] = t;
				tn.vec[1] = t;
				tn.vec[2] = t;
				planes= (co+remove_neg(dir))*(dims/ndims)-(dims/2);
				//for (int i = 0; i < NEWTON_STEPS; i++){
				for (int i=0; i < NEWTON_STEPS; i++){

					tn -= (a*tn*tn*tn+b*tn*tn+c*tn+d-planes)/(3*a*tn*tn+2*b*tn+c);
				}
				for (int i=0; i < 3; i++){
					if (abs(3*a.vec[i]*t*t+2*b.vec[i]*t+c.vec[i]) < 1e-8) tn.vec[i]=t+1;
				}
				//int min_index = argmin_not_nan(tn);
				int min_index = argmin(tn);
				//if (isnan(tn.vec[min_index]) ) return; // Things are really really really bad.
				if (t > tn.vec[min_index])tn.vec[min_index] = 0.0/0.0;
				t = tn.vec[min_index];



				p = a*t*t*t+b*t*t+c*t+d;
				co.vec[min_index] += dir.vec[min_index];

				id=co_to_idx(co,ndims);
				if (min(co) >= 0 && co < ndims ) atomicAdd(&(image[id]),norm(p_old-p)*proj);
				//if (min(co) >= 0 && co < ndims) image[id]=max((REAL)j,image[id]);
				//if (min(co) >= 0 && co < ndims) image[id]+=proj*length;
				p_old = p;


			}

	}

}
*/

template <class REAL> __global__ void Gadgetron::crop_splines_kernel(vector_td<REAL,3> * splines, REAL* projections, const  vector_td<REAL,3>  dims, const  vector_td<REAL,3>  origin, const int proj_dim,const REAL background,int offset)
{
	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x+offset;

	if (idx < proj_dim){
		const int sid = idx*4;

		const vector_td<REAL,3> half_dims = dims/((REAL)2.0);



		REAL t,told;

		//Load in points to registers
		vector_td<REAL,3> p0 = splines[sid]-origin;
		vector_td<REAL,3> p1 = splines[sid+1]-origin;
		vector_td<REAL,3> m0 = splines[sid+2];
		vector_td<REAL,3> m1 = splines[sid+3];

		REAL length = norm(p1-p0);
		m0 *= length/norm(m0);
		m1 *= length/norm(m1);
		vector_td<REAL,3> p,pt0,pt1;

		t=0;
		for (int i = 0; i < MAXSTEP; i++){
			told = t;
			t = ((REAL) i)/MAXSTEP;
			//t2 = t*t;

			//p = (2*t3-3*t2+1)*p0+(t3-2*t2+t)*m0+(3*t2-2*t3)*p1+(t3-t2)*m1+half_dims;
			p=t*m0+p0+half_dims;

			if ( min(p) >= 0 && p < dims) break;

		}
		t = told;
		//t2 = t*t;

		//pt0 =  (2*t3-3*t2+1)*p0+(t3-2*t2+t)*m0+(3*t2-2*t3)*p1+(t3-t2)*m1; //Calculate new starting point
		pt0=t*m0+p0;

		t = 0;
		for (int i = 0; i < MAXSTEP; i++){
				told = t;
				t = ((REAL) i)/MAXSTEP;
				//t2 = t*t;

				//p = (2*t3-3*t2+1)*p0+(t3-2*t2+t)*m0+(3*t2-2*t3)*p1+(t3-t2)*m1+half_dims;
				p=p1-t*m1+half_dims;
				if ( min(p) >= 0 && p < dims) break;

		}
		t = told;

		pt1=p1-t*m1;
		REAL deltaLength = norm(p1-pt1)+norm(p0-pt0);
		projections[idx] -= deltaLength*background;
		splines[sid]=pt0;
		splines[sid+1]=pt1;


	}
}


template <class REAL> __global__ void Gadgetron::rescale_directions_kernel(vector_td<REAL,3> * splines, REAL* projections, const  vector_td<REAL,3>  dims,  const int proj_dim, const int offset )
{
	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x+offset;

	if (idx < proj_dim){
		const int sid = idx*4;

		//Load in points to registers
		vector_td<REAL,3> p0 = splines[sid];
		vector_td<REAL,3> p1 = splines[sid+1];
		vector_td<REAL,3> m0 = splines[sid+2];
		vector_td<REAL,3> m1 = splines[sid+3];

		m0 /= norm(m0);
		m1 /= norm(m1);
		REAL length = norm(p1-p0);


		splines[sid+2]=m0*length;
		splines[sid+3]=m1*length;

	}
}

template <class REAL> __global__ void Gadgetron::points_to_coefficients(vector_td<REAL,3> * splines, int dim,int offset)
{


	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x+offset;
	const int sid = 4*idx;
	if (idx < dim){

		//Load in points to registers
		vector_td<REAL,3> p0 = splines[sid]; //Position at entrance
		vector_td<REAL,3> p1 = splines[sid+1]; // Position at exit
		vector_td<REAL,3> m0 = splines[sid+2]; // Direction at entrance
		vector_td<REAL,3> m1 = splines[sid+3]; // Direction at exit


		vector_td<REAL,3> a = 2*p0+m0+m1-2*p1;
		vector_td<REAL,3> b = -3*p0-2*m0+3*p1-m1;
		vector_td<REAL,3> c = m0;
		vector_td<REAL,3> d = p0;
		splines[sid]=a;
		splines[sid+1]=b;
		splines[sid+2]=c;
		splines[sid+3]=d;


	}

}

template <class REAL> __global__ void Gadgetron::spline_trapz_kernel(vector_td<REAL,3> * splines, REAL* lengths, int dim, int offset)
{
	//Admiral Ackbarz
	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x+offset;
	const int sid = 4*idx;
	if (idx < dim){
		REAL res = 0;
		REAL s1;
		//Load in points to registers
		vector_td<REAL,3> a = splines[sid];
		vector_td<REAL,3> b = splines[sid+1];
		vector_td<REAL,3> c = splines[sid+2];
		vector_td<REAL,3> d = splines[sid+3];

		REAL t = 0;
		REAL s0 = norm(d);

		for (int i = 1; i < INT_STEPS; i++){
			t = ((REAL) i)/INT_STEPS;
			s1 = norm(c+t*(2*b+t*3*a));
			res += (s0+s1)/(2*INT_STEPS);
			s0 = s1;
		}
		lengths[idx]=res;
	}

}

/***
 * The Hansen correctional facility for young cubic splines.
 * Corrects the path length with the ratio between the length of the straight line approximation and the length of the spline
 */
template <class REAL> __global__ void Gadgetron::length_correction_kernel(vector_td<REAL,3> * splines, REAL* projections, int dim, int offset)
{
	//Admiral Ackbarz
	const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x+offset;
	const int sid = 4*idx;
	if (idx < dim){
		REAL res = 0;
		REAL s1;
		//Load in points to registers
		vector_td<REAL,3> a = splines[sid];
		vector_td<REAL,3> b = splines[sid+1];
		vector_td<REAL,3> c = splines[sid+2];
		vector_td<REAL,3> d = splines[sid+3];

		REAL t = 0;

		REAL s0 = norm(d);

		for (int i = 1; i < INT_STEPS; i++){
			t = ((REAL) i)/INT_STEPS;
			s1 = norm(c+t*(2*b+t*3*a));
			res += (s0+s1)/(2*INT_STEPS);
			s0 = s1;
		}
		REAL tmp = norm(a+b+c)/res;
		projections[idx] *= tmp;
		if (idx-offset == 0)	printf("Correction: %f\n",tmp);
	}

}


template __global__ void Gadgetron::forward_kernel<float>(const float * __restrict__, float* ,const vector_td<float,3>  * __restrict__ ,  const vector_td<float,3> ,
		const typename intd<3>::Type, const int , const int );

template __global__ void Gadgetron::backwards_kernel<float>(const float* __restrict__ projections, float* __restrict__ image,
		const vector_td<float,3> * __restrict__ splines,  const vector_td<float,3> dims,
		const typename intd<3>::Type ndims, const int proj_dim, const int offset);

template __global__ void Gadgetron::crop_splines_kernel<float>(vector_td<float,3> * splines, float* projections, const  vector_td<float,3>  dims, const  vector_td<float,3>  origin,const int proj_dim,float background,int offset);
template __global__ void Gadgetron::rescale_directions_kernel<float>(vector_td<float,3> * splines, float* projections, const  vector_td<float,3>  dims,  const int proj_dim, const int offset);


template __global__ void Gadgetron::points_to_coefficients<float>(vector_td<float,3> * splines, int dim,int offset);

template __global__  void Gadgetron::length_correction_kernel<float>(vector_td<float,3> * splines, float* projections, int dim, int offset);
/*
template<> __global__ void forward_kernel<float>(float* image, float* projections,
		vector_td<float,3> * splines,  const vector_td<float,3> dims,
		const typename uintd<3>::Type ndims, const int proj_dim, const int offset);
		*/
/*


template<> __global__ void backwards_kernel<float>;
template<> __global__ void rescale_splines_kernel<float>;

template<> __global__ void points_to_coefficients<float>;

template<> __global__ void spline_trapz_kernel<float>;
*/
