
/*

-------------------------------------------------------------------------------------
Implementation details header. For internal use only in the cuda NFFT implementation.
-------------------------------------------------------------------------------------

*/


#ifndef _PREPROCESS_PRIVATE_HPP_
#define _PREPROCESS_PRIVATE_HPP_

#include <vector>

// Use circular kernel (as opposed to squared kernel)
//#define CIRCULAR_KERNEL

namespace mr_recon
{

	/*
		Domain priority queue for the NFFT_H preprocessing
	*/

	class pqueue_elem
	{
	public:

		pqueue_elem( unsigned int num_strips, unsigned int domain_index ){ 
			number_of_strips = num_strips;
			this->domain_index = domain_index;
		}

		unsigned int number_of_strips;
		unsigned int domain_index;		
	};

	class pqueue_compare 
	{ 
	public: 
		int operator()( const pqueue_elem &x, const pqueue_elem &y ) { 
			return x.number_of_strips < y.number_of_strips; 
		} 
	};


	/*
		A sample strip denotes a set of successive trajectory positions.
		It is represented by the index to the first sample and the number of samples.
	*/

	class sample_strip
	{
	public:

		sample_strip( unsigned int idx );

		unsigned int startIndex;
		unsigned int count;
	};


	/*

		A gridDomain denotes a sub-block of the d-dimensional Cartesian grid.
		Block dimensions are given in a RealIntArray representing two diaginal corners of the box.

		In 2D: Pt1 = (x1, y1), Pt2 = (x2, y2) : 		| x1 x2 |
														| y1 y2 |

		In 3D: Pt1 = (x1, y1, z1), Pt2 = (x2, y2, z2):	| x1 x2 |
														| y1 y2 |
														| z1 z2 |


		A gridDomain also contains a vector of sample strips.

	*/


	class gridDomain{
	public:

		gridDomain();
		gridDomain( mr_recon::RealFloatArray &dimensions );

		void set( mr_recon::RealFloatArray &dimensions );
		mr_recon::RealFloatArray& get();

		bool inside( double *pt );
		bool inside( float *pt );
		bool inside( unsigned int *pt );

#ifdef CIRCULAR_KERNEL
		void set( unsigned int *origin, float squared_radius );
		bool inside_circle( float *pt );
		unsigned int *origin;
		float squared_radius;
#endif
		mr_recon::RealFloatArray dimensions;
		std::vector<sample_strip> sample_strips; // Sample strip convolving into each domain
	};

}

#endif
