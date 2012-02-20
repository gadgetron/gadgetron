#pragma once

template< class T, unsigned int D > class vector_td
{
	public:
		T vec[D];
};

// Some typedefs for convenience (templated typedefs are not (yet) available in C++)

template< class REAL, unsigned int D > struct reald{
  typedef vector_td< REAL, D > Type;
};

template< unsigned int D > struct intd{
  typedef vector_td< int, D > Type;
};

template< unsigned int D > struct uintd{
  typedef vector_td< unsigned int, D > Type;
};

template< unsigned int D > struct floatd{
  typedef typename reald< float, D >::Type Type;
};

template< unsigned int D > struct doubled{
  typedef typename reald< double, D >::Type Type;
};



// Inherited structs with convenient constructors

struct intd1 : intd<1>::Type{
  typedef intd<1>::Type Type;
  intd1(){}
  intd1( int x ) { vec[0] = x; }
};

struct intd2 : intd<2>::Type{
  typedef intd<2>::Type Type;
  intd2(){}
  intd2( int x, int y ) { vec[0] = x; vec[1] = y; }
};

struct intd3 : intd<3>::Type{
  typedef intd<3>::Type Type;
  intd3(){}
  intd3( int x, int y, int z ) { vec[0] = x; vec[1] = y; vec[2] = z; }
};

struct intd4 : intd<4>::Type{
  typedef intd<4>::Type Type;
  intd4(){}
  intd4( int x, int y, int z, int w ) { vec[0] = x; vec[1] = y; vec[2] = z; vec[3] = w; }
};

struct uintd1 : uintd<1>::Type{
  typedef uintd<1>::Type Type;
  uintd1(){}
  uintd1( unsigned int x ) { vec[0] = x; }
};

struct uintd2 : uintd<2>::Type{
  typedef uintd<2>::Type Type;
  uintd2(){}
  uintd2( unsigned int x, unsigned int y ) { vec[0] = x; vec[1] = y; }
};

struct uintd3 : uintd<3>::Type{
  typedef uintd<3>::Type Type;
  uintd3(){}
  uintd3( unsigned int x, unsigned int y, unsigned int z ) { vec[0] = x; vec[1] = y; vec[2] = z; }
};

struct uintd4 : uintd<4>::Type{
  typedef uintd<4>::Type Type;
  uintd4(){}
  uintd4( unsigned int x, unsigned int y, unsigned int z, unsigned int w ) { vec[0] = x; vec[1] = y; vec[2] = z; vec[3] = w; }
};

struct floatd1 : floatd<1>::Type{
  typedef floatd<1>::Type Type;
  floatd1(){}
  floatd1( float x ) { vec[0] = x; }
};

struct floatd2 : floatd<2>::Type{
  typedef floatd<2>::Type Type;
  floatd2(){}
  floatd2( float x, float y ) { vec[0] = x; vec[1] = y; }
};

struct floatd3 : floatd<3>::Type{
  typedef floatd<3>::Type Type;
  floatd3(){}
  floatd3( float x, float y, float z ) { vec[0] = x; vec[1] = y; vec[2] = z; }
};

struct floatd4 : floatd<4>::Type{
  typedef floatd<4>::Type Type;
  floatd4(){}
  floatd4( float x, float y, float z, float w ) { vec[0] = x; vec[1] = y; vec[2] = z; vec[3] = w; }
};

struct doubled1 : doubled<1>::Type{
  typedef doubled<1>::Type Type;
  doubled1(){}
  doubled1( double x ) { vec[0] = x; }
};

struct doubled2 : doubled<2>::Type{
  typedef doubled<2>::Type Type;
  doubled2(){}
  doubled2( double x, double y ) { vec[0] = x; vec[1] = y; }
};

struct doubled3 : doubled<3>::Type{
  typedef doubled<3>::Type Type;
  doubled3(){}
  doubled3( double x, double y, double z ) { vec[0] = x; vec[1] = y; vec[2] = z; }
};

struct doubled4 : doubled<4>::Type{
  typedef doubled<4>::Type Type;
  doubled4(){}
  doubled4( double x, double y, double z, double w ) { vec[0] = x; vec[1] = y; vec[2] = z; vec[3] = w; }
};


