//
// Kaiser-Bessel convolution kernels
//

// 'bessi0' is taken from numerical recipes in C

__inline__ __device__ double 
bessi0(double x)
{
  double ax,ans,y;
  if ((ax=fabs(x))<3.75) 
    {
      y=x/3.75;
      y*=y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.0360768+y*0.0045813)))));
    } 
  else 
    {
      y=3.75/ax;
      ans=(-0.02057706+y*(0.02635537+y*(-0.01647633+(y*0.00392377))));
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.01328592+y*(0.00225319+y*(-0.00157565+y*(0.00916281+y*ans)))));
    }
  return ans;
}

__inline__ __device__ float 
bessi0(float x)
{
  float ax,ans,y;
  if ((ax=fabsf(x)) <3.75f) 
    {
      y=x/3.75f;
      y*=y;
      ans=1.0f+y*(3.5156229f+y*(3.0899424f+y*(1.2067492f+y*(0.2659732f+y*(0.0360768f+y*0.0045813f)))));
    } 
  else 
    {
      y=3.75f/ax;
      ans=(-0.02057706f+y*(0.02635537f+y*(-0.01647633f+(y*0.00392377f))));
      ans=(expf(ax)/sqrtf(ax))*(0.39894228f+y*(0.01328592f+y*(0.00225319f+y*(-0.00157565f+y*(0.00916281f+y*ans)))));
    }
  return ans;
}

// Kaiser Bessel according to Beatty et. al. IEEE TMI 2005;24(6):799-808.
// There is a slight difference wrt Jackson's formulation, IEEE TMI 1991;10(3):473-478.

__inline__ __device__ double
KaiserBessel( double u, double matrix_size_os, double one_over_W, double beta )
{
  double _tmp = 2.0*u*one_over_W;
  double tmp = _tmp*_tmp;
  double arg = beta*sqrt(1.0-tmp);
  double bessi = bessi0(arg);
  double ret = matrix_size_os*bessi*one_over_W;
  return ret;
}

__inline__ __device__ float
KaiserBessel( float u, float matrix_size_os, float one_over_W, float beta )
{
  float _tmp = 2.0f*u*one_over_W;
  float tmp = _tmp*_tmp;
  float arg = beta*sqrtf(1.0f-tmp);
  float bessi = bessi0(arg);
  float ret = matrix_size_os*bessi*one_over_W;
  return ret;
}

//
// Below the intended interface
//

template<class REALd, class REAL> __inline__ __device__ REAL
KaiserBessel( REALd u, REALd matrix_size_os, REAL one_over_W, REAL beta, uint2 fixedDims )
{
  REAL one; get_one(one);
  REAL phi_x = (fixedDims.x) ? one : KaiserBessel( u.x, matrix_size_os.x, one_over_W, beta );
  REAL phi_y = (fixedDims.y) ? one : KaiserBessel( u.y, matrix_size_os.y, one_over_W, beta );

  return phi_x*phi_y;
}

template<class REALd, class REAL> __inline__ __device__ REAL
KaiserBessel( REALd u, REALd matrix_size_os, REAL one_over_W, REAL beta, uint3 fixedDims )
{
  REAL one; get_one(one);
  REAL phi_x = (fixedDims.x) ? one : KaiserBessel( u.x, matrix_size_os.x, one_over_W, beta );
  REAL phi_y = (fixedDims.y) ? one : KaiserBessel( u.y, matrix_size_os.y, one_over_W, beta );
  REAL phi_z = (fixedDims.z) ? one : KaiserBessel( u.z, matrix_size_os.z, one_over_W, beta );

  return phi_x*phi_y*phi_z;
}

template<class REALd, class REAL> __inline__ __device__ REAL
KaiserBessel( REALd u, REALd matrix_size_os, REAL one_over_W, REAL beta, uint4 fixedDims )
{
  REAL one; get_one(one);
  REAL phi_x = (fixedDims.x) ? one : KaiserBessel( u.x, matrix_size_os.x, one_over_W, beta );
  REAL phi_y = (fixedDims.y) ? one : KaiserBessel( u.y, matrix_size_os.y, one_over_W, beta );
  REAL phi_z = (fixedDims.z) ? one : KaiserBessel( u.z, matrix_size_os.z, one_over_W, beta );
  REAL phi_w = (fixedDims.w) ? one : KaiserBessel( u.w, matrix_size_os.w, one_over_W, beta );

  return phi_x*phi_y*phi_z*phi_w;
}
