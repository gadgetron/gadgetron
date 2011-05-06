#ifndef _VEC_H
#define _VEC_H

#include "export.h"

namespace mr_recon_vectors
{
	template <class T> class DLLEXPORT vec2
	{
	public:
		vec2();
		vec2( const vec2<T>& );
		vec2( T, T );

		//	void normalize();
		//	T length() const;

		void set( T, T );
		void set( const vec2<T>& );

		T operator[]( unsigned int ) const;
		void operator= ( const vec2<T>& );
		void operator= ( int );

		vec2<T> operator+( const vec2<T>& ) const;
		vec2<T> operator-( const vec2<T>& ) const;
		vec2<T> operator*( T ) const;
		vec2<T> operator*( const vec2<T>& ) const;
		vec2<T> operator/( T ) const;
		vec2<T> operator/( const vec2<T>& ) const;

		void operator+=( const vec2<T>& );
		void operator-=( const vec2<T>& );
		void operator*=( T );
		void operator*=( const vec2<T>& );
		void operator/=( T );
		void operator/=( const vec2<T>& );

		bool operator== ( const vec2<T>& );

		T vec[2];
	};

	template <class T> class DLLEXPORT vec3
	{
	public:
		vec3();
		vec3( const vec3<T>& );
		vec3( T, T, T );

		//	void normalize();
		//	T length() const;

		void set( T, T, T );
		void set( const vec3<T>& );

		T operator[]( unsigned int ) const;
		void operator= ( const vec3<T>& );
		void operator= ( int );

		vec3<T> operator+( const vec3<T>& ) const;
		vec3<T> operator-( const vec3<T>& ) const;
		vec3<T> operator*( T ) const;
		vec3<T> operator*( const vec3<T>& ) const;
		vec3<T> operator/( T ) const;
		vec3<T> operator/( const vec3<T>& ) const;

		void operator+=( const vec3<T>& );
		void operator-=( const vec3<T>& );
		void operator*=( T );
		void operator*=( const vec3<T>& );
		void operator/=( T );
		void operator/=( const vec3<T>& );

		bool operator== ( const vec3<T>& );

		T vec[3];
	};

	template <class T> class DLLEXPORT vec4
	{
	public:
		vec4();
		vec4( const vec4<T>& );
		vec4( T, T, T, T );

		//	void normalize();
		//	T length() const;

		void set( T, T, T, T );
		void set( const vec4<T>& );

		T operator[]( unsigned int ) const;
		void operator= ( const vec4<T>& );
		void operator= ( int );

		vec4<T> operator+( const vec4<T>& ) const;
		vec4<T> operator-( const vec4<T>& ) const;
		vec4<T> operator*( T ) const;
		vec4<T> operator*( const vec4<T>& ) const;
		vec4<T> operator/( T ) const;
		vec4<T> operator/( const vec4<T>& ) const;

		void operator+=( const vec4<T>& );
		void operator-=( const vec4<T>& );
		void operator*=( T );
		void operator*=( const vec4<T>& );
		void operator/=( T );
		void operator/=( const vec4<T>& );

		bool operator== ( const vec4<T>& );

		T vec[4];
	};

	template <class T> class DLLEXPORT N_vec2 : public vec2<T>
	{
	public:
		N_vec2();
		N_vec2( const vec2<T>& );
		N_vec2( T, T );

		void operator= ( int );

		vec2<T> operator<<( const vec2<T>& ) const;
		vec2<T> operator>>( const vec2<T>& ) const;
		vec2<T> operator<<( T ) const;
		vec2<T> operator>>( T ) const;
		vec2<T> operator%( const vec2<T>& ) const;
		vec2<T> operator%( T ) const;
		void operator<<=( const vec2<T>& );
		void operator>>=( const vec2<T>& );
		void operator<<=( T );
		void operator>>=( T );
		void operator%=( const vec2<T>& );
		void operator%=( T );
		T prod() const;
	};

	template <class T> class DLLEXPORT N_vec3 : public vec3<T>
	{
	public:
		N_vec3();
		N_vec3( const vec3<T>& );
		N_vec3( T, T, T );

		void operator= ( int );

		vec3<T> operator<<( const vec3<T>& ) const;
		vec3<T> operator>>( const vec3<T>& ) const;
		vec3<T> operator<<( T ) const;
		vec3<T> operator>>( T ) const;
		vec3<T> operator%( const vec3<T>& ) const;
		vec3<T> operator%( T ) const;
		void operator<<=( const vec3<T>& );
		void operator>>=( const vec3<T>& );
		void operator<<=( T );
		void operator>>=( T );
		void operator%=( const vec3<T>& );
		void operator%=( T );
		T prod() const;
	};

	template <class T> class DLLEXPORT N_vec4 : public vec4<T>
	{
	public:
		N_vec4();
		N_vec4( const vec4<T>& );
		N_vec4( T, T, T, T );

		void operator= ( int );

		vec4<T> operator<<( const vec4<T>& ) const;
		vec4<T> operator>>( const vec4<T>& ) const;
		vec4<T> operator<<( T ) const;
		vec4<T> operator>>( T ) const;
		vec4<T> operator%( const vec4<T>& ) const;
		vec4<T> operator%( T ) const;
		void operator<<=( const vec4<T>& );
		void operator>>=( const vec4<T>& );
		void operator<<=( T );
		void operator>>=( T );
		void operator%=( const vec4<T>& );
		void operator%=( T );
		T prod() const;
	};
}

#endif
