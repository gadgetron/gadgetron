#include "vec.hpp"

using namespace mr_recon_vectors;

// VEC2

template<class T> 
vec2<T>::vec2()
{
}

template<class T> 
vec2<T>::vec2( const vec2<T>& from )
{
	set( from );
}

template<class T> 
vec2<T>::vec2( T x, T y )
{
	set( x, y );
}

/*
template<class T> void
vec2<T>::normalize()
{
    T len = sqrt( vec[0]*vec[0] + vec[1]*vec[1] );
	vec[0] /= len;
	vec[1] /= len;	
}
 
template<class T> T
vec2<T>::length() const
{
	return sqrt( vec[0]*vec[0] + vec[1]*vec[1] );
}
*/

template<class T> void
vec2<T>::set( T x, T y )
{
	vec[0] = x;
	vec[1] = y;	
}

template<class T> void
vec2<T>::set( const vec2<T>& from )
{
	vec[0] = from[0];
	vec[1] = from[1];
}
 
template<class T> T
vec2<T>::operator[]( unsigned int index ) const
{
	return vec[index];
}

template<class T> void
vec2<T>::operator= ( const vec2<T> &from )
{
	set( from );
}
	 
template<class T> void
vec2<T>::operator= ( int from )
{
	set( (T)from, (T)from );
}

template<class T> vec2<T>
vec2<T>::operator+( const vec2<T> &other ) const
{
	return vec2<T>( vec[0]+other[0], vec[1]+other[1] );
}

template<class T> vec2<T>
vec2<T>::operator-( const vec2<T> &other ) const
{
	return vec2<T>( vec[0]-other[0], vec[1]-other[1] );
}

template<class T> vec2<T>
vec2<T>::operator*( const vec2<T> &other ) const
{
	return vec2<T>( vec[0]*other[0], vec[1]*other[1] );
}

template<class T> vec2<T>
vec2<T>::operator/( const vec2<T> &other ) const
{
	return vec2<T>( vec[0]/other[0], vec[1]/other[1] );
}

template<class T> vec2<T>
vec2<T>::operator*( const T other ) const
{
	return vec2<T>( vec[0]*other, vec[1]*other );
}

template<class T> vec2<T>
vec2<T>::operator/( const T other ) const
{
	return vec2<T>( vec[0]/other, vec[1]/other );
}

template<class T> void
vec2<T>::operator+=( const vec2<T> &other )
{
	vec[0] += other[0];
	vec[1] += other[1];
}

template<class T> void
vec2<T>::operator-=( const vec2<T> &other )
{
	vec[0] -= other[0];
	vec[1] -= other[1];
}

template<class T> void
vec2<T>::operator*=( const vec2<T> &other )
{
	vec[0] *= other[0];
	vec[1] *= other[1];
}

template<class T> void
vec2<T>::operator/=( const vec2<T> &other )
{
	vec[0] /= other[0];
	vec[1] /= other[1];
}

template<class T> void
vec2<T>::operator*=( T other )
{
	vec[0] *= other;
	vec[1] *= other;
}

template<class T> void
vec2<T>::operator/=( T other )
{
	vec[0] /= other;
	vec[1] /= other;
}

template<class T> bool
vec2<T>::operator== ( const vec2<T> &from )
{
	return ((vec[0]==from[0])&&(vec[1]==from[1]));
}

// VEC3

template<class T> 
vec3<T>::vec3()
{
}

template<class T> 
vec3<T>::vec3( const vec3<T>& from )
{
	set( from );
}

template<class T> 
vec3<T>::vec3( T x, T y, T z )
{
	set( x, y, z );
}

/*
template<class T> void
vec3<T>::normalize()
{
    T len = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
	vec[0] /= len;
	vec[1] /= len;	
	vec[2] /= len;	
}
 
template<class T> T
vec3<T>::length() const
{
	return sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
}
*/
	
template<class T> void
vec3<T>::set( T x, T y, T z )
{
	vec[0] = x;
	vec[1] = y;	
	vec[2] = z;	
}

template<class T> void
vec3<T>::set( const vec3<T>& from )
{
	vec[0] = from[0];
	vec[1] = from[1];
	vec[2] = from[2];
}
 
template<class T> T
vec3<T>::operator[]( unsigned int index ) const
{
	return vec[index];
}

template<class T> void
vec3<T>::operator= ( const vec3<T> &from )
{
	set( from );
}
	 
template<class T> void
vec3<T>::operator= ( int from )
{
	set( (T)from, (T)from, (T)from );
}

template<class T> vec3<T>
vec3<T>::operator+( const vec3<T> &other ) const
{
	return vec3( vec[0]+other[0], vec[1]+other[1], vec[2]+other[2] );
}

template<class T> vec3<T>
vec3<T>::operator-( const vec3<T> &other ) const
{
	return vec3( vec[0]-other[0], vec[1]-other[1], vec[2]-other[2] );
}

template<class T> vec3<T>
vec3<T>::operator*( const vec3<T> &other ) const
{
	return vec3( vec[0]*other[0], vec[1]*other[1], vec[2]*other[2] );
}

template<class T> vec3<T>
vec3<T>::operator/( const vec3<T> &other ) const
{
	return vec3( vec[0]/other[0], vec[1]/other[1], vec[2]/other[2] );
}

template<class T> vec3<T>
vec3<T>::operator*( const T other ) const
{
	return vec3( vec[0]*other, vec[1]*other, vec[2]*other );
}

template<class T> vec3<T>
vec3<T>::operator/( const T other ) const
{
	return vec3( vec[0]/other, vec[1]/other, vec[2]/other );
}

template<class T> void
vec3<T>::operator+=( const vec3<T> &other )
{
	vec[0] += other[0];
	vec[1] += other[1];
	vec[2] += other[2];
}

template<class T> void
vec3<T>::operator-=( const vec3<T> &other )
{
	vec[0] -= other[0];
	vec[1] -= other[1];
	vec[2] -= other[2];
}

template<class T> void
vec3<T>::operator*=( const vec3<T> &other )
{
	vec[0] *= other[0];
	vec[1] *= other[1];
	vec[2] *= other[2];
}

template<class T> void
vec3<T>::operator/=( const vec3<T> &other )
{
	vec[0] /= other[0];
	vec[1] /= other[1];
	vec[2] /= other[2];
}

template<class T> void
vec3<T>::operator*=( T other )
{
	vec[0] *= other;
	vec[1] *= other;
	vec[2] *= other;
}

template<class T> void
vec3<T>::operator/=( T other )
{
	vec[0] /= other;
	vec[1] /= other;
	vec[2] /= other;
}

template<class T> bool
vec3<T>::operator== ( const vec3<T> &from )
{
	return ((vec[0]==from[0])&&(vec[1]==from[1])&&(vec[2]==from[2]));
}

// VEC4

template<class T> 
vec4<T>::vec4()
{
}

template<class T> 
vec4<T>::vec4( const vec4<T>& from )
{
	set( from );
}

template<class T> 
vec4<T>::vec4( T x, T y, T z, T w )
{
	set( x, y, z, w );
}

/*
template<class T> void
vec4<T>::normalize()
{
    T len = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] + vec[3]*vec[3] );
	vec[0] /= len;
	vec[1] /= len;	
	vec[2] /= len;	
	vec[3] /= len;	
}
 
template<class T> T
vec4<T>::length() const
{
	return sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] + vec[3]*vec[3] );
}
*/

template<class T> void
vec4<T>::set( T x, T y, T z, T w )
{
	vec[0] = x;
	vec[1] = y;	
	vec[2] = z;	
	vec[3] = w;	
}

template<class T> void
vec4<T>::set( const vec4<T>& from )
{
	vec[0] = from[0];
	vec[1] = from[1];
	vec[2] = from[2];
	vec[3] = from[3];
}
 
template<class T> T
vec4<T>::operator[]( unsigned int index ) const
{
	return vec[index];
}

template<class T> void
vec4<T>::operator= ( const vec4<T> &from )
{
	set( from );
}
	 
template<class T> void
vec4<T>::operator= ( int from )
{
	set( (T)from, (T)from, (T)from, (T)from );
}

template<class T> vec4<T>
vec4<T>::operator+( const vec4<T> &other ) const
{
	return vec4( vec[0]+other[0], vec[1]+other[1], vec[2]+other[2], vec[3]+other[3] );
}

template<class T> vec4<T>
vec4<T>::operator-( const vec4<T> &other ) const
{
	return vec4( vec[0]-other[0], vec[1]-other[1], vec[2]-other[2], vec[3]-other[3] );
}

template<class T> vec4<T>
vec4<T>::operator*( const vec4<T> &other ) const
{
	return vec4( vec[0]*other[0], vec[1]*other[1], vec[2]*other[2], vec[3]*other[3] );
}

template<class T> vec4<T>
vec4<T>::operator/( const vec4<T> &other ) const
{
	return vec4( vec[0]/other[0], vec[1]/other[1], vec[2]/other[2], vec[3]/other[3] );
}

template<class T> vec4<T>
vec4<T>::operator*( const T other ) const
{
	return vec4( vec[0]*other, vec[1]*other, vec[2]*other, vec[3]*other );
}

template<class T> vec4<T>
vec4<T>::operator/( const T other ) const
{
	return vec4( vec[0]/other, vec[1]/other, vec[2]/other, vec[3]/other );
}

template<class T> void
vec4<T>::operator+=( const vec4<T> &other )
{
	vec[0] += other[0];
	vec[1] += other[1];
	vec[2] += other[2];
	vec[3] += other[3];
}

template<class T> void
vec4<T>::operator-=( const vec4<T> &other )
{
	vec[0] -= other[0];
	vec[1] -= other[1];
	vec[2] -= other[2];
	vec[3] -= other[3];
}

template<class T> void
vec4<T>::operator*=( const vec4<T> &other )
{
	vec[0] *= other[0];
	vec[1] *= other[1];
	vec[2] *= other[2];
	vec[3] *= other[3];
}

template<class T> void
vec4<T>::operator/=( const vec4<T> &other )
{
	vec[0] /= other[0];
	vec[1] /= other[1];
	vec[2] /= other[2];
	vec[3] /= other[3];
}

template<class T> void
vec4<T>::operator*=( T other )
{
	vec[0] *= other;
	vec[1] *= other;
	vec[2] *= other;
	vec[3] *= other;
}

template<class T> void
vec4<T>::operator/=( T other )
{
	vec[0] /= other;
	vec[1] /= other;
	vec[2] /= other;
	vec[3] /= other;
}

template<class T> bool
vec4<T>::operator== ( const vec4<T> &from )
{
	return ((vec[0]==from[0])&&(vec[1]==from[1])&&(vec[2]==from[2])&&(vec[3]==from[3]));
}

// N_VEC2


template<class T> 
N_vec2<T>::N_vec2() : vec2<T>()
{
}

template<class T> 
N_vec2<T>::N_vec2( const vec2<T>& from ) : vec2<T>(from)
{
}

template<class T> 
N_vec2<T>::N_vec2( T x, T y ) : vec2<T>(x,y)
{
}

template<class T> void
N_vec2<T>::operator= ( int from )
{
	set( (T)from, (T)from );
}

template<class T> vec2<T>
N_vec2<T>::operator<<( const vec2<T>& v ) const
{
  	return N_vec2( this->vec[0]<<v[0], this->vec[1]<<v[1] );
}

template<class T> vec2<T>
N_vec2<T>::operator>>( const vec2<T>& v ) const
{
	return N_vec2( this->vec[0]>>v[0], this->vec[1]>>v[1] );
}

template<class T> vec2<T>
N_vec2<T>::operator<<( T v ) const
{
	return N_vec2( this->vec[0]<<v, this->vec[1]<<v );
}

template<class T> vec2<T>
N_vec2<T>::operator>>( T v ) const
{
	return N_vec2( this->vec[0]>>v, this->vec[1]>>v );
}

template<class T> void 
N_vec2<T>::operator<<=( const vec2<T>& v )
{
	this->vec[0]<<=v[0];
	this->vec[1]<<=v[1];
}

template<class T> void 
N_vec2<T>::operator>>=( const vec2<T>& v )
{
	this->vec[0]>>=v[0];
	this->vec[1]>>=v[1];
}

template<class T> void 
N_vec2<T>::operator<<=( T v )
{
	this->vec[0]<<=v;
	this->vec[1]<<=v;
}

template<class T> void 
N_vec2<T>::operator>>=( T v )
{
	this->vec[0]>>=v;
	this->vec[1]>>=v;
}

template<class T> vec2<T> 
N_vec2<T>::operator%( const vec2<T>& v ) const
{
	return N_vec2( this->vec[0]%v[0], this->vec[1]%v[1] );
}

template<class T> vec2<T> 
N_vec2<T>::operator%( T v ) const
{
	return vec2<T>( this->vec[0]%v, this->vec[1]%v );
}

template<class T> void
N_vec2<T>::operator%=( const vec2<T>& v )
{
	this->vec[0]%=v[0];
	this->vec[1]%=v[1];
}

template<class T> void 
N_vec2<T>::operator%=( T v )
{
	this->vec[0]%=v;
	this->vec[1]%=v;
}
template<class T> T 
N_vec2<T>::prod() const
{
	return this->vec[0]*this->vec[1];
}


// N_VEC3


template<class T> 
N_vec3<T>::N_vec3() : vec3<T>()
{
}

template<class T> 
N_vec3<T>::N_vec3( const vec3<T>& from ) : vec3<T>(from)
{
}

template<class T> 
N_vec3<T>::N_vec3( T x, T y, T z ) : vec3<T>(x,y,z)
{
}

template<class T> void
N_vec3<T>::operator= ( int from )
{
	set( (T)from, (T)from, (T)from );
}

template<class T> vec3<T>
N_vec3<T>::operator<<( const vec3<T>& v ) const
{
	return N_vec3( this->vec[0]<<v[0], this->vec[1]<<v[1], this->vec[2]<<v[2] );
}

template<class T> vec3<T>
N_vec3<T>::operator>>( const vec3<T>& v ) const
{
	return N_vec3( this->vec[0]>>v[0], this->vec[1]>>v[1], this->vec[2]>>v[2] );
}

template<class T> vec3<T>
N_vec3<T>::operator<<( T v ) const
{
	return N_vec3( this->vec[0]<<v, this->vec[1]<<v, this->vec[2]<<v );
}

template<class T> vec3<T>
N_vec3<T>::operator>>( T v ) const
{
	return N_vec3( this->vec[0]>>v, this->vec[1]>>v, this->vec[2]>>v );
}

template<class T> void 
N_vec3<T>::operator<<=( const vec3<T>& v )
{
	this->vec[0]<<=v[0];
	this->vec[1]<<=v[1];
	this->vec[2]<<=v[2];
}

template<class T> void 
N_vec3<T>::operator>>=( const vec3<T>& v )
{
	this->vec[0]>>=v[0];
	this->vec[1]>>=v[1];
	this->vec[2]>>=v[2];
}

template<class T> void 
N_vec3<T>::operator<<=( T v )
{
	this->vec[0]<<=v;
	this->vec[1]<<=v;
	this->vec[2]<<=v;
}

template<class T> void 
N_vec3<T>::operator>>=( T v )
{
	this->vec[0]>>=v;
	this->vec[1]>>=v;
	this->vec[2]>>=v;
}

template<class T> vec3<T> 
N_vec3<T>::operator%( const vec3<T>& v ) const
{
	return N_vec3( this->vec[0]%v[0], this->vec[1]%v[1], this->vec[2]%v[2] );
}

template<class T> vec3<T> 
N_vec3<T>::operator%( T v ) const
{
	return N_vec3( this->vec[0]%v, this->vec[1]%v, this->vec[2]%v );
}

template<class T> void
N_vec3<T>::operator%=( const vec3<T>& v )
{
	this->vec[0]%=v[0];
	this->vec[1]%=v[1];
	this->vec[2]%=v[2];
}

template<class T> void 
N_vec3<T>::operator%=( T v )
{
	this->vec[0]%=v;
	this->vec[1]%=v;
	this->vec[2]%=v;
}

template<class T> T 
N_vec3<T>::prod() const
{
	return this->vec[0]*this->vec[1]*this->vec[2];
}


// N_vec4


template<class T> 
N_vec4<T>::N_vec4() : vec4<T>()
{
}

template<class T> 
N_vec4<T>::N_vec4( const vec4<T>& from ) : vec4<T>(from)
{
}

template<class T> 
N_vec4<T>::N_vec4( T x, T y, T z, T w ) : vec4<T>(x,y,z,w)
{
}

template<class T> void
N_vec4<T>::operator= ( int from )
{
	set( (T)from, (T)from, (T)from, (T)from );
}

template<class T> vec4<T>
N_vec4<T>::operator<<( const vec4<T>& v ) const
{
	return N_vec4( this->vec[0]<<v[0], this->vec[1]<<v[1], this->vec[2]<<v[2], this->vec[3]<<v[3] );
}

template<class T> vec4<T>
N_vec4<T>::operator>>( const vec4<T>& v ) const
{
	return N_vec4( this->vec[0]>>v[0], this->vec[1]>>v[1], this->vec[2]>>v[2], this->vec[3]>>v[3] );
}

template<class T> vec4<T>
N_vec4<T>::operator<<( T v ) const
{
	return N_vec4( this->vec[0]<<v, this->vec[1]<<v, this->vec[2]<<v, this->vec[3]<<v );
}

template<class T> vec4<T>
N_vec4<T>::operator>>( T v ) const
{
	return N_vec4( this->vec[0]>>v, this->vec[1]>>v, this->vec[2]>>v, this->vec[3]>>v );
}

template<class T> void 
N_vec4<T>::operator<<=( const vec4<T>& v )
{
	this->vec[0]<<=v[0];
	this->vec[1]<<=v[1];
	this->vec[2]<<=v[2];
	this->vec[3]<<=v[3];
}

template<class T> void 
N_vec4<T>::operator>>=( const vec4<T>& v )
{
	this->vec[0]>>=v[0];
	this->vec[1]>>=v[1];
	this->vec[2]>>=v[2];
	this->vec[3]>>=v[3];
}

template<class T> void 
N_vec4<T>::operator<<=( T v )
{
	this->vec[0]<<=v;
	this->vec[1]<<=v;
	this->vec[2]<<=v;
	this->vec[3]<<=v;
}

template<class T> void 
N_vec4<T>::operator>>=( T v )
{
	this->vec[0]>>=v;
	this->vec[1]>>=v;
	this->vec[2]>>=v;
	this->vec[3]>>=v;
}

template<class T> vec4<T> 
N_vec4<T>::operator%( const vec4<T>& v ) const
{
	return N_vec4( this->vec[0]%v[0], this->vec[1]%v[1], this->vec[2]%v[2], this->vec[3]%v[3] );
}

template<class T> vec4<T> 
N_vec4<T>::operator%( T v ) const
{
	return N_vec4( this->vec[0]%v, this->vec[1]%v, this->vec[2]%v, this->vec[3]%v );
}

template<class T> void
N_vec4<T>::operator%=( const vec4<T>& v )
{
	this->vec[0]%=v[0];
	this->vec[1]%=v[1];
	this->vec[2]%=v[2];
	this->vec[3]%=v[3];
}

template<class T> void 
N_vec4<T>::operator%=( T v )
{
	this->vec[0]%=v;
	this->vec[1]%=v;
	this->vec[2]%=v;
	this->vec[3]%=v;
}

template<class T> T 
N_vec4<T>::prod() const
{
	return this->vec[0]*this->vec[1]*this->vec[2]*this->vec[3];
}


/* The declarations below are necessary to make sure that specific functions are included in the library */
template class DLLEXPORT vec2< double >;
template class DLLEXPORT vec3< double >;
template class DLLEXPORT vec4< double >;
template class DLLEXPORT vec2< float >;
template class DLLEXPORT vec3< float >;
template class DLLEXPORT vec4< float >;
template class DLLEXPORT vec2< int >;
template class DLLEXPORT vec3< int >;
template class DLLEXPORT vec4< int >;
template class DLLEXPORT vec2< unsigned int >;
template class DLLEXPORT vec3< unsigned int >;
template class DLLEXPORT vec4< unsigned int >;
template class DLLEXPORT vec2< char >;
template class DLLEXPORT vec3< char >;
template class DLLEXPORT vec4< char >;
template class DLLEXPORT N_vec2< int >;
template class DLLEXPORT N_vec3< int >;
template class DLLEXPORT N_vec4< int >;
template class DLLEXPORT N_vec2< unsigned int >;
template class DLLEXPORT N_vec3< unsigned int >;
template class DLLEXPORT N_vec4< unsigned int >;
template class DLLEXPORT N_vec2< char >;
template class DLLEXPORT N_vec3< char >;
template class DLLEXPORT N_vec4< char >;
