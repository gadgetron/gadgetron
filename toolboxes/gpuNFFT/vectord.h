#pragma once

template< class T, unsigned int D > struct vectord
{
  T vec[D];
};

//
// Integer types
//

template< unsigned int D > struct intd  : vectord< int, D >{};
template< unsigned int D > struct uintd : vectord< unsigned int, D >{};

struct __align__(8) intd2 : intd<2>{
  intd2(){}
  intd2( int x, int y ) { vec[0] = x; vec[1] = y; }
};

struct intd3 : intd<3>{
  intd3(){}
  intd3( int x, int y, int z ) { vec[0] = x; vec[1] = y; vec[2] = z; }
};

struct __align__(16) intd4 : intd<4>{
  intd4(){}
  intd4( int x, int y, int z, int w ) { vec[0] = x; vec[1] = y; vec[2] = z; vec[3] = w; }
};

struct __align__(8) uintd2 : uintd<2>{
  uintd2(){}
  uintd2( unsigned int x, unsigned int y ) { vec[0] = x; vec[1] = y; }
};

struct uintd3 : uintd<3>{
  uintd3(){}
  uintd3( unsigned int x, unsigned int y, unsigned int z ) { vec[0] = x; vec[1] = y; vec[2] = z; }
};

struct __align__(16) uintd4 : uintd<4>{
  uintd4(){}
  uintd4( unsigned int x, unsigned int y, unsigned int z, unsigned int w ) { vec[0] = x; vec[1] = y; vec[2] = z; vec[3] = w; }
};

//
// Real types
//

template< unsigned int D > struct floatd  : vectord< float, D > {};
template< unsigned int D > struct doubled : vectord< double, D > {};

struct __align__(8) floatd2 : floatd<2>{
  floatd2(){}
  floatd2( float x, float y ) { vec[0] = x; vec[1] = y; }
};

struct floatd3 : floatd<3>{
  floatd3(){}
  floatd3( float x, float y, float z ) { vec[0] = x; vec[1] = y; vec[2] = z; }
};

struct __align__(16) floatd4 : floatd<4>{
  floatd4(){}
  floatd4( float x, float y, float z, float w ) { vec[0] = x; vec[1] = y; vec[2] = z; vec[3] = w; }
};

struct __align__(16) doubled2 : doubled<2>{
  doubled2(){}
  doubled2( double x, double y ) { vec[0] = x; vec[1] = y; }
};

struct __align__(16) doubled3 : doubled<3>{
  doubled3(){}
  doubled3( double x, double y, double z ) { vec[0] = x; vec[1] = y; vec[2] = z; }
};

struct __align__(16) doubled4 : doubled<4>{
  doubled4(){}
  doubled4( double x, double y, double z, double w ) { vec[0] = x; vec[1] = y; vec[2] = z; vec[3] = w; }
};

//
// Real complex types
//

template< class REAL > struct real_complex : vectord< REAL, 2 > {};

struct __align__(8) float_complex : real_complex<float> {};
struct __align__(16) double_complex : real_complex<double> {};
