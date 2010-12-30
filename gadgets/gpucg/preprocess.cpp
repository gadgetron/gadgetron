//
// Generate the necessary data structures for the CUDA convolution in the NFFT and NFFT_H algorithms
//

#include "preprocess.hpp"
#include "preprocess_private.hpp"
#include "preprocess_private.hcu"
#include "NFFT.hcu"
#include "uint_util.hcu"
#include "float_util.hcu"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <queue>
#include <iostream>

#include <vector_types.h>
#include <vector_functions.h>

using namespace std;
using namespace mr_recon;

// MACRO definitions

// Convert 1D index 'i' to d-dimensional index 'res' corresponding to dSizes.
#define index_to_d_dim( i, d, dSizes, res )			\
{								\
	unsigned int _denominator_ = 1;				\
	for( unsigned int _d_=0; _d_<d; _d_++ ){		\
		res[_d_] = (i/_denominator_)%dSizes[_d_];	\
		_denominator_ *= dSizes[_d_];			\
	}							\
}

#define _BIG_FLOAT 99999999999999999999999.9f;


// Cartesian grid encoding direction (longest dimension for strips)
enum DIR { XDIR=0, YDIR, ZDIR, TDIR };

/*

	------------------------
	Template implementation. 
	------------------------

	All the work is left for the constructors of the individual classes.

*/

template< class UINTd, class FLOATd, char TYPE > NFFT_plan<UINTd, FLOATd, TYPE>*
mr_recon::preprocess_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, RealFloatArray *trajectory )
{

	NFFT_plan<UINTd, FLOATd, TYPE> *tmp = new NFFT_plan<UINTd, FLOATd, TYPE>( matrix_size, matrix_size_os, fixed_dims, domain_size_samples, domain_size_coils, W );

	tmp->preprocess( trajectory );

	if( tmp->successfully_preprocessed )
		return tmp;
	else{
		delete tmp;
		return 0x0;
	}
}

template< class UINTd, class FLOATd, char TYPE > NFFT_H_plan<UINTd, FLOATd, TYPE>*
mr_recon::preprocess_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_coils, float W, RealFloatArray *trajectory )
{
	NFFT_H_plan<UINTd, FLOATd, TYPE> *tmp =  new NFFT_H_plan<UINTd, FLOATd, TYPE>( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_coils, W );

	tmp->preprocess( trajectory );

	if( tmp->successfully_preprocessed )
		return tmp;
	else{
		delete tmp;
		return 0x0;
	}
}

template< class UINTd, class FLOATd, char TYPE > NFFT_iteration_plan<UINTd, FLOATd, TYPE>*
mr_recon::preprocess_NFFT( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, RealFloatArray *trajectory )
{
	NFFT_iteration_plan<UINTd, FLOATd, TYPE> *tmp =  new NFFT_iteration_plan<UINTd, FLOATd, TYPE>( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_samples, domain_size_coils, W );

	tmp->preprocess( trajectory );

	if( tmp->successfully_preprocessed )
		return tmp;
	else{
		delete tmp;
		return 0x0;
	}
}

/*
	NFFT preprocessing implementation
	---------------------------------

	Our goal is to build the data structure in linear time O(M+N) where M is the number of samples and N is the number of grid cells.
	For this we need to utilize data structures with constant (O(1)) time access for each element.

	For each sample domain (a sequence of samples) we generate a vector of bit vectors to mark which grid cells convolve onto the domain.
		- the outmost vector size denote the multiplum of the (d-1) "narrowest" grid dimensions.
		- the (inner) bit vector the marks whether each individual grid cell in the "longest" grid dimension contributes to the convolution.
		- each sequential series of grid points in the long dimension will be considered a strip, encoded by a dD coordinate, a direction vector, and the number of samples.

*/

template< class UINTd, class FLOATd, char TYPE > 
NFFT_plan< UINTd, FLOATd, TYPE >::NFFT_plan( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, unsigned int domain_size_samples, unsigned int domain_size_coils, float W )
{
	this->matrix_size = matrix_size;
	this->matrix_size_os = matrix_size_os;

	this->fixed_dims = fixed_dims;

	this->stripsDir_NFFT = 0x0;
	this->stripOrigins_NFFT = 0x0;

	this->sample_positions = 0x0;

	this->stripsDirDevPtr_NFFT = 0x0;
	this->stripOriginsDevPtr_NFFT = 0x0;
	this->stripsMapDevPtr_NFFT = 0x0;
	this->stripLengthsDevPtr_NFFT = 0x0;

	this->stripsMap_NFFT = 0x0;

	this->domain_size_samples = domain_size_samples;
	this->domain_count_samples = 0;

	if( domain_size_coils>8 ) // Limit of shared memory
		this->domain_size_coils = 8;
	else
		this->domain_size_coils = domain_size_coils;

	this->number_of_coils = domain_size_coils;

	// Determine dimensionality 'd'
	this->d = sizeof(UINTd)/sizeof(unsigned int);
	assert( (sizeof(UINTd)%sizeof(unsigned int)) == 0 );
	assert(	d>1 && d<=4 );

	// If there is any fixed dimension then make sure there is no oversampling here. Also the corresponding grid_size entry will be set to one.
	this->alpha = 1.0f; 
	for( unsigned int i=0; i<d; i++ ){
		if( ((unsigned int*)&fixed_dims)[i] ){
			if( ((unsigned int*)&matrix_size)[i] != ((unsigned int*)&matrix_size_os)[i] ){
				printf("\nWarning: Oversampling factor MUST be one for all fixed dimensions. Enforcing this.\n");
				((unsigned int*)&(this->matrix_size_os))[i] = ((unsigned int*)&matrix_size)[i];
			}
		}
		else{
			// Determine oversampling factor 'alpha'
			this->alpha = (float)((unsigned int*)&matrix_size_os)[0]/(float)((unsigned int*)&matrix_size)[0];
		}
	}

	// Check matrix size consistency
	for( unsigned int _d=0; _d<d; _d++ ){
		if( !((unsigned int*)&fixed_dims)[_d] && !(fabs(alpha*((unsigned int*)&matrix_size)[_d] - (float)(((unsigned int*)&matrix_size_os)[_d]))<0.001f) ){
			printf("\nFATAL ERROR: Oversampling factor cannot vary between non-fixed dimensions. Quitting.\n");
			exit(1);
		}
	}

	this->W = W;

	if( W <= 0.0f ){
		printf("\nWarning: 'W' must be strictly positive. Resetting to 1.0f.\n");
		this->W = 1.0f;
	}

	this->number_of_strips_NFFT = 0;
	this->number_of_samples = 0;

	this->stripLengths_NFFT = 0x0;

	this->sample_positions_DevPtr = 0x0;
		
	this->samples_per_projection = 0;
	this->projections_per_frame = 0;
	this->angular_offset = 0;
	this->frames_per_rotation = 0;
	this->interframe_rotation = 0.0f;
	this->gc_factor = 0.0f;
	this->total_projections_f = 0.0f;

	this->max_strips_per_domain_NFFT = 0;

	this->successfully_preprocessed = false;
	this->initialized = false;
}

template< class UINTd, class FLOATd, char TYPE > bool 
NFFT_plan< UINTd, FLOATd, TYPE >::preprocess( RealFloatArray *trajectory )
{	
	
	successfully_preprocessed = false;

	// Set number of samples and number of domains
	number_of_samples = 1;
	for( int i=0; i<trajectory->get_number_of_dimensions()-1; i++ ){
	  number_of_samples *= trajectory->get_size(i);
	}

	domain_count_samples = number_of_samples/domain_size_samples;

	if( !((trajectory->get_number_of_dimensions()==2)||(trajectory->get_number_of_dimensions()==3)) ){
		printf("\nERROR: Number of dimensions in trajectory array is not two or three.\n");
		return false;
	}

	if( (trajectory->get_number_of_dimensions()==2) && ((unsigned int)trajectory->get_size(1)!=d) ){
		printf("\nERROR: Height of trajectory array is not 'd'.\n");
		return false;
	}

	if( (trajectory->get_number_of_dimensions()==3) && ((unsigned int)trajectory->get_size(2)!=d-1) ){
		printf("\nERROR: Height of trajectory array is not 'd'-1.\n");
		return false;
	}
	
	if( number_of_samples%domain_size_samples ) {
	        printf("\nERROR: We currently require the number of samples to be a multiple of the domain length.\n");
		return false;				
	}

	if( (number_of_samples/domain_size_samples)%2 ){
	        printf("\nERROR: We currently require the number of samples to be an EVEN multiple of the domain length.\n");
		return false;				
	}
	

	// Allocate array memory

	if( stripsMap_NFFT )
		delete[] stripsMap_NFFT;

	if( stripsDir_NFFT )
		delete[] stripsDir_NFFT;

	stripsMap_NFFT = new uint2[domain_count_samples];
	stripsDir_NFFT = new UINTd[domain_count_samples];

	if( !stripsMap_NFFT || !stripsDir_NFFT )
	{
		printf("\nERROR: Cannot allocate memory for utility arrays (NFFT)!\n");
		return false;	
	}

	// We operate on positve numbers [0;matrix_size_os] whereas trajectory positions are centered arounf the origin
	UINTd bias = matrix_size_os>>1;

	// Fix bias for fixed uneven dimensions
	for( unsigned int _d=0; _d<d; _d++ )
		if( ((unsigned int*)&fixed_dims)[_d] && (((unsigned int*)&matrix_size_os)[_d]%2) )
			((unsigned int*)&bias)[_d]++;

	// Kernel and half kernel width (oversampled!)
	const float aW = alpha * W;
	const float aW_half = aW / 2.0f;

	// Allocate memory for domain bitvectors
	vector<vector<bool> >* domains_bitvectors = new vector<vector<bool> >[domain_count_samples];

	// Offsets vectors (required to restore grid positions)
	vector<UINTd>* domains_offsets = new vector<UINTd>[domain_count_samples];

	// Directions vector for strips
	vector<DIR> dirs(domain_count_samples);

	// Make samples array

	if( sample_positions )
		delete[] sample_positions;

	sample_positions = new FLOATd[number_of_samples];

	if( !sample_positions ){
		printf("\nERROR: out of memory during preprocessing.\n");
		return false;
	}
	
	for( unsigned int sample=0; sample<number_of_samples; sample++ )
	{
		// d-dimensional sample position 'x' of NDArray template type T
		FLOATd x;

		for( unsigned int i=0; i<d; i++ ){

		  if( ((unsigned int*)&fixed_dims)[i] ){
		    if( trajectory->get_number_of_dimensions()==2)
		      ((float*)&x)[i] = (trajectory->get_data_ptr())[i*number_of_samples+sample] + ((unsigned int*)&bias)[i];
		    else
		      ((float*)&x)[i] = floorf((float)(((unsigned int*)&matrix_size)[i])*((float)(sample/trajectory->get_size(0)))/((float)trajectory->get_size(1)));
		  }
		  else
		    ((float*)&x)[i] = alpha*((trajectory->get_data_ptr())[i*number_of_samples+sample]) + ((unsigned int*)&bias)[i];
		  
		  // Make sure we have not exceeded our matrix size?
		  assert( (((float*)&x)[i]>=0) && (((float*)&x)[i]<=((unsigned int*)&matrix_size_os)[i]) );
		}
		
		// Update the sample_positions with point 'x'.
		sample_positions[sample] = x;
	}


	// Variables used in the loop below.

	FLOATd minx;						// Domain bounding box minumum
	FLOATd maxx;						// Domain bounding box maximum
	FLOATd deltax;						// Domain bounding box delta (maxx - minx)
	unsigned int *narrow_entries = new unsigned int[d-1];	// Number of entries in the narrow dimensions.
	unsigned int num_narrow_entries;			// Total number of entries in the narrow dimensions.
	unsigned int num_long_entries = 0;				// Total number of entries in the longest dimension.
	unsigned int *gridPos_offset = new unsigned int[d-1];	// Grid "minimum" coordinate in the narrow dimensions
	unsigned int gridPos_offset_long = 0;			// Grid "minimum" coordinate in the longest dimension
	unsigned int *gridPos_domain = new unsigned int[d-1];	// Grid position in the domain bounding box in the narrow dimensions

	// Per sample values
	FLOATd x;						// Sample position
	uint2 *gridInterval_narrow = new uint2[d-1];		// First and last grid coordinate in each narrow dimension
	uint2 gridInterval_long = make_uint2(0,0);				// First and last grid coordinate in longest dimension (sample encoding direction)
	unsigned int *gridDelta_narrow = new unsigned int[d-1];	// Number of grid points in each narrow dimension
	unsigned int gridDelta_long;				// Number of grid points in in longest dimension (sample encoding direction)
	unsigned int *gridPos_narrow = new unsigned int[d-1];	// Current grid coordinate narrow dimensions (for iteration)
	unsigned int gridPos_long;				// Current grid coordinate long dimension (for iteration)
	unsigned int *_offset = new unsigned int[d];

	// Process each sample domain, ie. 'domain_length_NFFT' samples, and find grid points that convolve onto either sample

	for( unsigned int i=0; i<domain_count_samples; i++ ){
	
		for( unsigned int _d=0; _d<d; _d++ ){
			((float*)&minx)[_d] = _BIG_FLOAT;
			((float*)&maxx)[_d] = 0.0f;
		}

		// Determine which direction produce the longest strips

		for( unsigned int j=0; j<domain_size_samples; j++ ){

			x = sample_positions[domain_size_samples*i+j];

			for( unsigned int _d=0; _d<d; _d++ ){

				if( ((float*)&x)[_d]<((float*)&minx)[_d] )
					((float*)&minx)[_d] = ((float*)&x)[_d];

				if( ((float*)&x)[_d]>((float*)&maxx)[_d] )
					((float*)&maxx)[_d] = ((float*)&x)[_d];
			}
		}

		for( unsigned int _d=0; _d<d; _d++ )
			((float*)&deltax)[_d] = ((float*)&maxx)[_d]-((float*)&minx)[_d];

		// Set 'dir' to the direction to the _largest_ span.
		
		DIR dir = (DIR)0; // Overwritten below, but set initially to avoid compiler warnings
		float max_delta = -1.0f;

		for( unsigned int _d=0; _d<d; _d++ ){
			if( ((float*)&deltax)[_d]>max_delta ){
				max_delta = ((float*)&deltax)[_d];
				dir = (DIR) _d;
			}
		}
		
		// We need 'dir' later to construct grid strips
		dirs[i] = dir;
		
		if( ((unsigned int*)&fixed_dims)[dir] ){
			printf("\nPreprocessing implementation error!!! %.1f, %d\n", max_delta, dir);
			return false;
		}

		// Multiply the narrow dimensions to determine the outer vector length
		
		num_narrow_entries = 1;

		for( unsigned int _d=0, _add = 0; _d<d; _d++ ){

			int offset = (int) ceil( ((float*)&minx)[_d]-aW_half );
			int _end = (int) floor( ((float*)&maxx)[_d]+aW_half );
			unsigned int length = (unsigned int) _end-offset+1;

			if( (unsigned int)dir != _d ){
				if( ((unsigned int*)&fixed_dims)[_d] )
					gridPos_offset[_d-_add] = (unsigned int) floor( ((float*)&minx)[_d] );
				else
					gridPos_offset[_d-_add] = (unsigned int) max( 0, offset );

				if( gridPos_offset[_d-_add]+length > ((unsigned int*)&matrix_size_os)[_d] )
					length = ((unsigned int*)&matrix_size_os)[_d]-gridPos_offset[_d-_add];

				if( ((unsigned int*)&fixed_dims)[_d] )
					narrow_entries[_d-_add] = 1;
				else
					narrow_entries[_d-_add] = (offset<0) ? length+offset : length;

				num_narrow_entries *= narrow_entries[_d-_add];
			}
			else{
				gridPos_offset_long = (unsigned int) max( 0, offset );
				if( gridPos_offset_long+length > ((unsigned int*)&matrix_size_os)[dir] )
					length = ((unsigned int*)&matrix_size_os)[dir]-gridPos_offset_long;
				num_long_entries = (offset<0) ? length+offset : length;
				_add++;
			}
		}

		// Resize vectors for domain 'i'.

		domains_bitvectors[i].resize( num_narrow_entries );
		domains_offsets[i].resize( num_narrow_entries );

		for( unsigned int entry=0; entry<num_narrow_entries; entry++ ){
			(domains_bitvectors[i])[entry].resize( num_long_entries, false );

			UINTd offset;
			index_to_d_dim( entry, d-1, narrow_entries, _offset )

			for( unsigned int _d=0, _add=0; _d<d; _d++ ){
			  if( _d != (unsigned int)dir )
					((unsigned int*)&offset)[_d] = gridPos_offset[_d-_add]+_offset[_d-_add];
				else{
					((unsigned int*)&offset)[_d] = gridPos_offset_long;
					_add++;
				}
			}

			(domains_offsets[i])[entry] = offset;
		}


		// Now we are ready to fill the domain vectros

		for( unsigned int j=0; j<domain_size_samples; j++ ){
		
			// Get sample position
		  for( unsigned int _d=0; _d<d; _d++ )
		    if( ((unsigned int*)&fixed_dims)[_d] ){
		      
		      if( trajectory->get_number_of_dimensions()==2)		 
			((float*)&x)[_d] = trajectory->get_data_ptr()[_d*number_of_samples+domain_size_samples*i+j] + ((unsigned int*)&bias)[_d];
		      else
			((float*)&x)[_d] = floorf((float)(((unsigned int*)&matrix_size)[_d])*((float)((domain_size_samples*i+j)/trajectory->get_size(0)))/((float)trajectory->get_size(1)));
		    }
		    else
		      ((float*)&x)[_d] = alpha*(trajectory->get_data_ptr()[_d*number_of_samples+domain_size_samples*i+j]) + ((unsigned int*)&bias)[_d];
		  
		  // Determine which grid points to iterate over
		  unsigned int grid_pts_count_narrow = 1, grid_pts_count_long = 0;
		  for( unsigned int _d=0, _add=0; _d<d; _d++ ){
			  if( _d != (unsigned int)dir ){
					// TODO: are we recomputing what was already done above?
					if( ((unsigned int*)&fixed_dims)[_d] ){
						gridInterval_narrow[_d-_add].x = gridInterval_narrow[_d-_add].y = (unsigned int) floor( ((float*)&x)[_d] );
					}
					else{
						gridInterval_narrow[_d-_add].x = (unsigned int) max( 0, (int) ceil( ((float*)&x)[_d] - aW_half ));
						gridInterval_narrow[_d-_add].y = min( ((unsigned int*)&matrix_size_os)[_d]-1, (unsigned int) floor( ((float*)&x)[_d] + aW_half ));
					}
					gridDelta_narrow[_d-_add] = gridInterval_narrow[_d-_add].y-gridInterval_narrow[_d-_add].x+1;
					grid_pts_count_narrow *= gridDelta_narrow[_d-_add];
				}
				else{
					gridInterval_long.x = (unsigned int) max( 0, (int) ceil( ((float*)&x)[dir] - aW_half ));
					gridInterval_long.y = min( ((unsigned int*)&matrix_size_os)[_d]-1, (unsigned int) floor( ((float*)&x)[dir] + aW_half ));
					gridDelta_long = gridInterval_long.y-gridInterval_long.x+1;
					grid_pts_count_long = gridDelta_long;
					_add++;
				}
			}

			// Iterate over kernel
			for( unsigned int k=0; k<grid_pts_count_narrow; k++ ){

				// We need the (d-1)D narrow dimensions coord
				unsigned int _mul = gridDelta_narrow[0];
				gridPos_narrow[0] = (k%_mul)+gridInterval_narrow[0].x;
				for( unsigned int _d=1; _d<d-1; _d++ ){
					gridPos_narrow[_d] = (k/_mul)+(gridInterval_narrow[_d]).x;
					_mul *= gridDelta_narrow[_d];
				}

#ifdef _DEBUG
				for( unsigned int _d=0, _add=0; _d<d; _d++ ){
					if( _d != dir ){
						assert( gridPos_narrow[_d-_add] >= 0 );
						assert( gridPos_narrow[_d-_add] < ((unsigned int*)&matrix_size_os)[_d] );
					}
					else
						_add++;
				}
#endif

				for( unsigned int l=0; l<grid_pts_count_long; l++ ){

					gridPos_long = gridInterval_long.x+l;

					// Index 'k' is local to the sample point. Find the corresponding domain-vector index ('gridPos_domain_idx').
					unsigned int gridPos_domain_idx = 0;
					for( unsigned int _d=0; _d<d-1; _d++ ){
						unsigned int diff = gridPos_narrow[_d]-gridPos_offset[_d];
						assert(diff<narrow_entries[_d]); 
						gridPos_domain[_d] = diff;	
					}
					unsigned int _mul = 1;
					for( unsigned int _d=0; _d<d-1; _d++ ){
						gridPos_domain_idx += gridPos_domain[_d]*_mul;
						_mul *= narrow_entries[_d];
					}

					assert( (gridPos_long-gridPos_offset_long)<num_long_entries );
					((domains_bitvectors[i])[gridPos_domain_idx])[gridPos_long-gridPos_offset_long] = true;
				}
			}
		}
	}

	vector<UINTd> _stripOrigins;
	vector<unsigned int> _stripLengths;

	unsigned int accNumStrips = 0;
	bool inside = false;

	for( unsigned int i=0; i<domain_count_samples; i++ ){

		unsigned int accNumStrips_save = accNumStrips;

		for( unsigned int j=0; j<domains_bitvectors[i].size(); j++ ){
			for( unsigned int k=0; k<(domains_bitvectors[i])[j].size(); k++ ){
				if( ((domains_bitvectors[i])[j])[k] ){
					if( !inside ){
						UINTd origin = (domains_offsets[i])[j];
						((unsigned int*) &origin)[dirs[i]] += k;
						_stripOrigins.push_back(origin);
						_stripLengths.push_back(1);
						accNumStrips++;
						inside = true;
					}
					else{
						_stripLengths[accNumStrips-1]++;
					}
				}
				else
					inside = false;
			}
			inside = false;
		}
	
		stripsMap_NFFT[i] = make_uint2( accNumStrips_save, accNumStrips-accNumStrips_save );
		uint4 noCompilerWarning = make_uint4( dirs[i]==XDIR, dirs[i]==YDIR, dirs[i]==ZDIR, dirs[i]==TDIR );
		uint4* noCompilerWarningPtr = &noCompilerWarning;
		stripsDir_NFFT[i] = *((UINTd*)noCompilerWarningPtr);
	}

	number_of_strips_NFFT = accNumStrips;
	
	if( stripOrigins_NFFT )
		delete[] stripOrigins_NFFT;

	if( stripLengths_NFFT )
		delete[] stripLengths_NFFT;

	stripOrigins_NFFT = new UINTd[number_of_strips_NFFT];
	stripLengths_NFFT = new unsigned int[number_of_strips_NFFT];

	if( !stripOrigins_NFFT || !stripLengths_NFFT )
	{
		printf("\nERROR: Cannot allocate memory for utility arrays (NFFT)!\n");
		return false;	
	}

	assert( (_stripOrigins.size() == _stripLengths.size()) && (_stripOrigins.size() == number_of_strips_NFFT) );
	memcpy( stripOrigins_NFFT, &_stripOrigins[0], number_of_strips_NFFT*sizeof(UINTd) );
	memcpy( stripLengths_NFFT, &_stripLengths[0], number_of_strips_NFFT*sizeof(unsigned int) );


	// Clean up
	delete[] domains_bitvectors;
	delete[] domains_offsets;
	delete[] narrow_entries;
	delete[] gridPos_offset;
	delete[] gridPos_domain;
	delete[] gridInterval_narrow;
	delete[] gridDelta_narrow;
	delete[] gridPos_narrow;
	delete[] _offset;

	successfully_preprocessed = true;

	return true;
}

/*
	NFFT_H preprocessing implementation
*/

template< class UINTd, class FLOATd, char TYPE > 
NFFT_H_plan< UINTd, FLOATd, TYPE >::NFFT_H_plan( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_coils, float W )
{
	this->matrix_size = matrix_size;
	this->matrix_size_os = matrix_size_os;

	this->domain_size_grid = domain_size_grid;

	this->fixed_dims = fixed_dims;

	sample_positions = 0x0;

	this->stripsMap_NFFT_H = 0x0;
	this->strips_NFFT_H = 0x0;
	this->domainsMapDevPtr_NFFT_H = 0x0;
	this->stripsMapDevPtr_NFFT_H = 0x0;
	this->stripsDevPtr_NFFT_H = 0x0;

	// Determine dimensionality 'd'
	this->d = sizeof(UINTd)/sizeof(unsigned int);
	assert( (sizeof(UINTd)%sizeof(unsigned int)) == 0 );
	assert(	d>1 && d<=4 );

	// If there is any fixed dimension then make sure there is no oversampling here. Also the corresponding grid_size entry will be set to one.
	this->alpha = 1.0f; 
	for( unsigned int i=0; i<d; i++ ){
		if( ((unsigned int*)&fixed_dims)[i] ){
			if( ((unsigned int*)&matrix_size)[i] != ((unsigned int*)&matrix_size_os)[i] ){
				printf("\nWarning: Oversampling factor MUST be one for all fixed dimensions. Enforcing this.\n");
				((unsigned int*)&(this->matrix_size_os))[i] = ((unsigned int*)&matrix_size)[i];
			}
			if( ((unsigned int*)&domain_size_grid)[i] != 1 ){
				printf( "\nWarning: Domain_size_grid MUST be one for all fixed dimensions. Enforcing this.\n" );
				((unsigned int*)&(this->domain_size_grid))[i] = 1;
			}
		}
		else{
			// Determine oversampling factor 'alpha'
			this->alpha = (float)((unsigned int*)&matrix_size_os)[0]/(float)((unsigned int*)&matrix_size)[0];
		}
	}

	// Check matrix size consistency
	for( unsigned int _d=0; _d<d; _d++ ){
		if( !((unsigned int*)&fixed_dims)[_d] && !(fabs(alpha*((unsigned int*)&matrix_size)[_d] - (float)(((unsigned int*)&matrix_size_os)[_d]))<0.001f) ){
			printf("\nFATAL ERROR: Oversampling factor cannot vary between non-fixed dimensions. Quitting.\n");
			exit(1);
		}
	}

	this->W = W;

	if( W <= 0.0f ){
		printf("\nWarning: 'W' must be positive. Resetting 'W' to 1.0f.\n");
		this->W = 1.0f;
	}

	
	UINTd zeros; uint_to_uintd( 0, 0, &zeros );
	UINTd ones; uint_to_uintd( 1, 1, &ones );

	// Determine neccessary wrap size for "folding" convolution. Make sure it is a multiple of the domain_size. No wrapping in fixed dimensions.
	if( weak_eq( domain_size_grid, zeros )){
		printf("\nERROR: domain_size_grid cannot be 0 in any dimension. Quitting.\n");
		exit(1);
	}
	matrix_size_wrap = (((((unsigned int) ceilf(this->W/2.0f))<<1)+(domain_size_grid-ones))/domain_size_grid)*domain_size_grid;

	for( unsigned int _d=0; _d<d; _d++ ){
		if( ((unsigned int*)&fixed_dims)[_d] && ((unsigned int*)&matrix_size_wrap)[_d] ){
			((unsigned int*)&matrix_size_wrap)[_d] = 0;
		}
	}

	this->domain_count_grid = (matrix_size_os+matrix_size_wrap)/domain_size_grid;

	if( domain_size_coils>8 ) // Limit of shared memory
		this->domain_size_coils = 8;
	else
		this->domain_size_coils = domain_size_coils;

	this->number_of_coils = domain_size_coils;

	this->calculate_deapodization = true;
	this->do_deapodization = true;
	this->clean_deapodization = true;

	this->domainsMap_NFFT_H = 0x0;

	this->number_of_strips_NFFT_H = 0;
	this->number_of_samples = 0;

	this->sample_positions_DevPtr = 0x0;
	this->weights_DevPtr = 0x0;

	this->samples_per_projection = 0;
	this->projections_per_frame = 0;
	this->angular_offset = 0;
	this->frames_per_rotation = 0;
	this->interframe_rotation = 0.0f;
	this->gc_factor = 0.0f;
	this->total_projections_f = 0.0f;

	for( unsigned int i=0; i<d; i++ ){

		if( ((unsigned int*)&fixed_dims)[i] )
			continue;

		float _domain_count = (float)(((unsigned int*)(&this->matrix_size_os))[i]+((unsigned int*)(&this->matrix_size_wrap))[i])/(float)((unsigned int*)(&this->domain_size_grid))[i];
		if( _domain_count != ((unsigned int) _domain_count) ){
			printf("\nFATAL ERROR: We currently require the oversampled grid size to be a multiple of the domain size. Quitting.\n");
			exit(1);				
		}
		if( ((unsigned int) _domain_count)%2 ){
			printf("\nFATAL ERROR: We currently require the oversampled grid size to be an EVEN multiple of the domain size. Quitting.\n");
			exit(1);				
		}
	}

	this->traj_positions = 0x0;
	this->tuples_last = 0x0;
	this->bucket_begin = this->bucket_end = 0x0;

	this->successfully_preprocessed = false;
	this->initialized = false;
}

template< class UINTd, class FLOATd, char TYPE > bool
NFFT_H_plan< UINTd, FLOATd, TYPE >::preprocess( RealFloatArray *trajectory )
{
	
	successfully_preprocessed = false;

	// Set number of samples
	number_of_samples = 1;
	for( int i=0; i<trajectory->get_number_of_dimensions()-1; i++ ){
	  number_of_samples *= trajectory->get_size(i);
	}

	// Check input validity

	if( !((trajectory->get_number_of_dimensions()==2)||(trajectory->get_number_of_dimensions()==3)) ){
		printf("\nERROR: Number of dimensions in trajectory array is not two or three.\n");
		return false;
	}

	if( (trajectory->get_number_of_dimensions()==2) && ((unsigned int)trajectory->get_size(1)!=d) ){
		printf("\nERROR: Height of trajectory array is not 'd'.\n");
		return false;
	}

	// We operate on positve numbers [0;matrix_size_os+matrix_size_wrap]
	UINTd bias = (matrix_size_os+matrix_size_wrap)>>1;

	// Fix bias for fixed uneven dimensions
	for( unsigned int _d=0; _d<d; _d++ ){
		if( ((unsigned int*)&fixed_dims)[_d] && (((unsigned int*)&matrix_size_os)[_d]%2) ){
			((unsigned int*)&bias)[_d]++;
			printf("\nWARNING: We never really tested uneven matrix sizes (in a fixed dimension)!!!\n");
		}
	}

	if( sample_positions )
		delete[] sample_positions;
	
	sample_positions = new FLOATd[number_of_samples];

	if( !sample_positions ){
		printf("\nERROR: out of memory during preprocessing.\n");
		return false;
	}

	unsigned int number_of_domains_NFFT_H = prod(domain_count_grid);

	// Allocate memory for domain arrays
	gridDomain *domains = new gridDomain[number_of_domains_NFFT_H]();

	if( !domains ){
		printf("\nERROR: out of memory allocating domain data!\n");
		return false;
	}

	// Initialize gridDomains (all samples within this box/cube will be processed by the corresponding thread)

	UINTd domIdx;
#ifdef CIRCULAR_KERNEL
	UINTd dom_origin;
#endif
	unsigned int *domIdxPtr = &domIdx.x;
	FLOATd farCorner;

	for( unsigned int i=0; i<d; i++ )
		((float*)&farCorner)[i] = -1.0f*((float)(((unsigned int*)&domain_size_grid)[i]*(((unsigned int*)&domain_count_grid)[i]>>1)));

	float margin = alpha*W/2.0f;

	// Calculate RealFloatArray providing boundaries for gridDomain 'i'.
	const unsigned int gridDomainDimSizes[2] = { 2, d };
	RealFloatArray dims( 2, (int*)gridDomainDimSizes );	
	float *data = dims.get_data_ptr();

	unsigned int *dummy_domain_count_grid = new unsigned int[d];
	for( unsigned int i=0; i<d; i++ )
		dummy_domain_count_grid[i] = ((unsigned int*)&domain_count_grid)[i];

	for( unsigned int i=0; i<number_of_domains_NFFT_H; i++ ){

		// Get domain i's d-dimensional index
		index_to_d_dim( i, d, dummy_domain_count_grid, domIdxPtr );

		for( unsigned int j=0; j<d; j++ ){
			if( ((unsigned int*)&fixed_dims)[j] ){
				data[2*j+0] = ((unsigned int*)&domIdx)[j]*((unsigned int*)&domain_size_grid)[j]-0.5f; 
				data[2*j+1] = ((unsigned int*)&domIdx)[j]*((unsigned int*)&domain_size_grid)[j]+0.5f;
			}
			else{
				data[2*j+0] = ((float*)&farCorner)[j]-margin+(float)(((unsigned int*)&domIdx)[j]*((unsigned int*)&domain_size_grid)[j]+((unsigned int*)&bias)[j]); 
				data[2*j+1] = ((float*)&farCorner)[j]+margin-1.0f+(float)((((unsigned int*)&domIdx)[j]+1)*((unsigned int*)&domain_size_grid)[j]+((unsigned int*)&bias)[j]); 
			}

#ifdef CIRCULAR_KERNEL
			dom_origin.vec[j] = domIdx[j]*domain_size_grid[j];
#endif
		}

		domains[i].set( dims );

#ifdef CIRCULAR_KERNEL
 		domains[i]->set( (unsigned int*) &dom_origin, margin );
#endif
	}


	/*
		Now determine which samples convolve into which gridDomains
	*/

	// 'corner' variables are used to obtain linear complexity by only looping over the domains our convolution window actually reaches
	FLOATd corner1, corner2;

	// Some variables used in the loop below
	unsigned int numSampleDomains = 1;
	UINTd sampleDomMin, sampleDomMax, sampleDomDelta;

	unsigned int *dummy_sampleDomDelta = new unsigned int[d];

	// Append samples to the grid domains they convolve into

	for( unsigned int sample=0; sample<(unsigned int) number_of_samples; sample++ )
	{
		// d-dimensional sample position 'x' of NDArray template type T
		FLOATd x;

		for( unsigned int i=0; i<d; i++ ){
		  if( ((unsigned int*)&fixed_dims)[i] ){
		    if( trajectory->get_number_of_dimensions()==2)
		      ((float*)&x)[i] = (trajectory->get_data_ptr())[i*number_of_samples+sample] + ((unsigned int*)&bias)[i];
		    else
		      ((float*)&x)[i] = floorf((float)(((unsigned int*)&matrix_size)[i])*((float)(sample/trajectory->get_size(0)))/((float)trajectory->get_size(1)));
		  }
		  else
		    ((float*)&x)[i] = alpha*((trajectory->get_data_ptr())[i*number_of_samples+sample]) + ((unsigned int*)&bias)[i];
		}
		
		// Update the sample_positions with point 'x'.
		sample_positions[sample] = x;

		// Update 'corner' util vectors before processing each sample
		for( unsigned int i=0; i<d; i++ ){

			((float*)&corner1)[i] = (((float*)&x)[i]-((float*)&farCorner)[i]-margin-((unsigned int*)&bias)[i])/((unsigned int*)&domain_size_grid)[i]-1;
			((float*)&corner2)[i] = (((float*)&x)[i]-((float*)&farCorner)[i]+margin-((unsigned int*)&bias)[i])/((unsigned int*)&domain_size_grid)[i];

			if( ((float*)&corner1)[i] < 0.0f ) 
				((float*)&corner1)[i] = 0.0f;
			if( ((float*)&corner2)[i] > (((unsigned int*)&domain_count_grid)[i]-1.0f) )
				((float*)&corner2)[i] = (((unsigned int*)&domain_count_grid)[i]-1.0f);
		}

		// Calculate the number of grid domains 'numSampleDomains' we must iterate over.
		// For a given setup there is a maximum value for 'numSampleDomains', i.e. linear complexity in the O() notation.

		numSampleDomains = 1;

		for( unsigned int i=0; i<d; i++ ){
			if( ((unsigned int*)&fixed_dims)[i] ){
				// Don't contribute in the domain search
			  ((unsigned int*)&sampleDomMin)[i] = ((unsigned int*)&sampleDomMax)[i] = (unsigned int)(((float*)&x)[i]/((unsigned int*)&domain_size_grid)[i]);
				((unsigned int*)&sampleDomDelta)[i] = 1;
			}
			else{
				((unsigned int*)&sampleDomMin)[i] = (unsigned int)((float*)&corner1)[i];
				((unsigned int*)&sampleDomMax)[i] = (unsigned int)((float*)&corner2)[i];
				((unsigned int*)&sampleDomDelta)[i] = ((unsigned int*)&sampleDomMax)[i]-((unsigned int*)&sampleDomMin)[i]+1;
				numSampleDomains *= ((unsigned int*)&sampleDomDelta)[i];
			}
		}

		// Loop over domains to determine which sample position 'x' is inside.

		for( unsigned int i=0; i<numSampleDomains; i++ ){

			for( unsigned int ii=0; ii<d; ii++ )
				dummy_sampleDomDelta[ii] = ((unsigned int*)&sampleDomDelta)[ii];

			index_to_d_dim( i, d, dummy_sampleDomDelta, domIdxPtr );

			// Find global domain index
			unsigned int globalIdx = 0, mul = 1;
			for( unsigned int j=0; j<d; j++ ){
				globalIdx += (((unsigned int*)&sampleDomMin)[j]+((unsigned int*)&domIdx)[j])*mul;
				mul *= ((unsigned int*)&domain_count_grid)[j];
			}

			gridDomain *domain = &domains[globalIdx];

			// Check if sample 'x' is inside the domain with index 'globalIdx'
			if( domain->inside((float*)&x) ){	// TODO: what if x was a double? FIXME!

				// When using circular kernel we need an additional check
#ifdef CIRCULAR_KERNEL
				if( !domain->inside_circle((float*)&x) )
					continue;
#endif
				// Each domain has sample_strip vector. If there are no strips at this point we are certainly about to start a new.
				size_t num_strips = domain->sample_strips.size();

				// Check if we need to make a new strip or append the current sample to an existing strip
				if( num_strips >0 && ( domain->sample_strips[num_strips-1].startIndex + domain->sample_strips[num_strips-1].count == sample) ){

					// Append to existing strip
					domain->sample_strips[num_strips-1].count += 1;
				}
				else{
					// Start new strip
					domain->sample_strips.push_back( sample_strip(sample) );
				}
			}
		}
	}


	// gridDomains have now been computed.
	// Extract vecctor data to build the class arrays.

	size_t totalNumStrips=0, maxNumStrips=0, maxSampleStartIdx = 0, maxSampleCount = 0;

	for( unsigned int i=0; i<number_of_domains_NFFT_H; i++ ){
		size_t numStrips = domains[i].sample_strips.size();
		totalNumStrips += numStrips;
		if( numStrips>maxNumStrips )
			maxNumStrips = numStrips;

		for( unsigned int strip=0; strip<numStrips; strip++ ){
			unsigned int sampleStartIdx = domains[i].sample_strips[strip].startIndex;
			unsigned int sampleCount = domains[i].sample_strips[strip].count;
			if(  sampleStartIdx > maxSampleStartIdx )
				maxSampleStartIdx = sampleStartIdx;
			if(  sampleCount > maxSampleCount )
				maxSampleCount = sampleCount;
		}
	}


	// Fill output object

	number_of_strips_NFFT_H = (unsigned int) totalNumStrips;
	
	if( domainsMap_NFFT_H )
		delete[] domainsMap_NFFT_H;

	if( stripsMap_NFFT_H )
		delete[] stripsMap_NFFT_H;

	if( strips_NFFT_H )
		delete[] strips_NFFT_H;

	domainsMap_NFFT_H = new unsigned int[number_of_domains_NFFT_H];
	stripsMap_NFFT_H =	new uint2[number_of_domains_NFFT_H];
	strips_NFFT_H = new uint2[number_of_strips_NFFT_H];

	if( !domainsMap_NFFT_H || !stripsMap_NFFT_H || !strips_NFFT_H )
	{
		printf("\nERROR: Cannot allocate memory for utility arrays!.\n");
		return false;	
	}


	// Fill arrays
	
	unsigned int accNumStrips = 0;
#ifndef _DEBUG
	priority_queue<pqueue_elem, vector<pqueue_elem>, pqueue_compare> pqueue;
#endif
	for( unsigned int i=0; i<number_of_domains_NFFT_H; i++ ){

		size_t numStrips = domains[i].sample_strips.size();

		((unsigned int*)&(stripsMap_NFFT_H[i]))[0] = (unsigned int) accNumStrips;
		((unsigned int*)&(stripsMap_NFFT_H[i]))[1] = (unsigned int) numStrips;
	
#ifndef _DEBUG
		pqueue.push( pqueue_elem((unsigned int)numStrips,i) );
#endif
		for( unsigned int strip=0; strip<numStrips; strip++ ){

			((unsigned int*)&(strips_NFFT_H[accNumStrips]))[0] = domains[i].sample_strips[strip].startIndex;
			((unsigned int*)&(strips_NFFT_H[accNumStrips]))[1] = domains[i].sample_strips[strip].count;

			accNumStrips++;
		}
	}

	// Build domainsMap from priority queue
#ifndef _DEBUG
	for( unsigned int i=0; i<number_of_domains_NFFT_H; i++ ){
		domainsMap_NFFT_H[i] = pqueue.top().domain_index;
		pqueue.pop();
	}
	assert( pqueue.empty() );
#endif


	// Clean up
	delete[] dummy_domain_count_grid;
	delete[] dummy_sampleDomDelta;
	delete[] domains;

	successfully_preprocessed = true;

	return true;
}

/*
	NFFT_iteration preprocessing implementation
*/

template< class UINTd, class FLOATd, char TYPE > 
NFFT_iteration_plan< UINTd, FLOATd, TYPE >::NFFT_iteration_plan( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_samples, unsigned int domain_size_coils, float W )
{
	// Constructor for iteration plan

	this->matrix_size = matrix_size;
	this->matrix_size_os = matrix_size_os;

	this->W = W;

	if( W <= 0.0f ){
		printf("\nWarning: 'W' must be strictly positive. Resetting to 1.0f.\n");
		this->W = 1.0f;
	}

	// Determine dimensionality 'd'
	this->d = sizeof(UINTd)/sizeof(unsigned int);
	assert( (sizeof(UINTd)%sizeof(unsigned int)) == 0 );
	assert(	d>1 && d<=4 );

	this->fixed_dims = fixed_dims;

	this->domain_size_grid = domain_size_grid;

	UINTd zeros; uint_to_uintd( 0, 0, &zeros );
	UINTd ones; uint_to_uintd( 1, 1, &ones );

	// Determine neccessary wrap size for "folding" convolution. Make sure it is a multiple of the domain_size. No wrapping in fixed dimensions.
	if( weak_eq( domain_size_grid, zeros )){
		printf("\nERROR: domain_size_grid cannot be 0 in any dimension. Quitting.\n");
		exit(1);
	}
	matrix_size_wrap = (((((unsigned int) ceilf(this->W/2.0f))<<1)+(domain_size_grid-ones))/domain_size_grid)*domain_size_grid;

	for( unsigned int _d=0; _d<d; _d++ ){
		if( ((unsigned int*)&fixed_dims)[_d] && ((unsigned int*)&matrix_size_wrap)[_d] ){
			((unsigned int*)&matrix_size_wrap)[_d] = 0;
		}
	}

	this->domain_count_grid = (matrix_size_os+matrix_size_wrap)/domain_size_grid;

	this->stripsDir_NFFT = 0x0;
	this->stripOrigins_NFFT = 0x0;

	this->sample_positions = 0x0;

	this->stripsDirDevPtr_NFFT = 0x0;
	this->stripOriginsDevPtr_NFFT = 0x0;
	this->stripsMapDevPtr_NFFT = 0x0;
	this->stripLengthsDevPtr_NFFT = 0x0;

	this->stripsMap_NFFT = 0x0;
	this->stripsMap_NFFT_H = 0x0;
	this->strips_NFFT_H = 0x0;

	this->domainsMapDevPtr_NFFT_H = 0x0;
	this->stripsMapDevPtr_NFFT_H = 0x0;
	this->stripsDevPtr_NFFT_H = 0x0;

	this->domain_size_samples = domain_size_samples;
	this->domain_count_samples = 0;

	if( domain_size_coils>8 ) // Limit of shared memory
		this->domain_size_coils = 8;
	else
		this->domain_size_coils = domain_size_coils;

	this->number_of_coils = domain_size_coils;

	// If there is any fixed dimension then make sure there is no oversampling here. Also the corresponding grid_size entry will be set to one.
	this->alpha = 1.0f; 
	for( unsigned int i=0; i<d; i++ ){
		if( ((unsigned int*)&fixed_dims)[i] ){
			if( ((unsigned int*)&matrix_size)[i] != ((unsigned int*)&matrix_size_os)[i] ){
				printf("\nWarning: Oversampling factor MUST be one for all fixed dimensions. Enforcing this.\n");
				((unsigned int*)&(this->matrix_size_os))[i] = ((unsigned int*)&matrix_size)[i];
			}
			if( ((unsigned int*)&domain_size_grid)[i] != 1 ){
				printf( "\nWarning: Domain_size_grid MUST be one for all fixed dimensions. Enforcing this.\n" );
				((unsigned int*)&(this->domain_size_grid))[i] = 1;
			}
		}
		else{
			// Determine oversampling factor 'alpha'
			this->alpha = (float)matrix_size_os.x/(float)matrix_size.x;
		}
	}

	// Check matrix size consistency
	for( unsigned int _d=0; _d<d; _d++ ){
		if( !((unsigned int*)&fixed_dims)[_d] && !(fabs(alpha*((unsigned int*)&(this->matrix_size))[_d] - (float)(((unsigned int*)&(this->matrix_size_os))[_d]))<0.001f) ){
			printf("\nFATAL ERROR: Oversampling factor cannot vary between non-fixed dimensions. Quitting.\n");
			exit(1);
		}
	}

	this->number_of_samples = 0;
	this->number_of_strips_NFFT = 0;
//	this->number_of_threads_NFFT_H = 0;
	this->number_of_strips_NFFT_H = 0;

	this->domainsMap_NFFT_H = 0x0;
	this->stripLengths_NFFT = 0x0;

	this->sample_positions_DevPtr = 0x0;
	this->weights_DevPtr = 0x0;

	this->samples_per_projection = 0;
	this->projections_per_frame = 0;
	this->angular_offset = 0;
	this->frames_per_rotation = 0;
	this->interframe_rotation = 0.0f;
	this->gc_factor = 0.0f;
	this->total_projections_f = 0.0f;

	this->max_strips_per_domain_NFFT = 0;

	this->regularizationDevPtr = 0x0;
	this->CSMdevPtr = 0x0;
	this->intensity_correction_magnitudes_image_DevPtr = 0x0;

	this->traj_positions = 0x0;
	this->tuples_last = 0x0;
	this->bucket_begin = this->bucket_end = 0x0;

	successfully_preprocessed = false;
	initialized = false;
}


template< class UINTd, class FLOATd, char TYPE > bool
NFFT_iteration_plan< UINTd, FLOATd, TYPE >::preprocess( RealFloatArray *trajectory )
{
	NFFT_plan< UINTd, FLOATd, TYPE > *pre_NFFT = new NFFT_plan< UINTd, FLOATd, TYPE >( matrix_size, matrix_size_os, fixed_dims, domain_size_samples, domain_size_coils, W );
	NFFT_H_plan< UINTd, FLOATd, TYPE > *pre_NFFT_H = new NFFT_H_plan< UINTd, FLOATd, TYPE >( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_coils, W );

	this->successfully_preprocessed = false;

	pre_NFFT->preprocess( trajectory );
	pre_NFFT_H->preprocess( trajectory );

	if( !pre_NFFT->successfully_preprocessed || !pre_NFFT->successfully_preprocessed ){
		NFFT_cleanup(&pre_NFFT);
		NFFT_cleanup(&pre_NFFT_H);
		return false;
	}

	assert( pre_NFFT->number_of_samples == pre_NFFT_H->number_of_samples );
	this->number_of_samples = pre_NFFT->number_of_samples;

	this->domain_count_samples = pre_NFFT->domain_count_samples;

	if( this->sample_positions )
		delete[] this->sample_positions;

	this->sample_positions = new FLOATd[number_of_samples];
	memcpy( sample_positions, pre_NFFT->sample_positions, number_of_samples*sizeof(FLOATd) );

	this->number_of_strips_NFFT = pre_NFFT->number_of_strips_NFFT;

	if( this->stripsMap_NFFT )
		delete[] this->stripsMap_NFFT;

	this->stripsMap_NFFT = new uint2[domain_count_samples];
	memcpy( stripsMap_NFFT, pre_NFFT->stripsMap_NFFT, domain_count_samples*sizeof(uint2) );

	if( this->stripsDir_NFFT )
		delete[] this->stripsDir_NFFT;

	this->stripsDir_NFFT = new UINTd[domain_count_samples];
	memcpy( stripsDir_NFFT, pre_NFFT->stripsDir_NFFT, domain_count_samples*sizeof(UINTd) );

	if( this->stripOrigins_NFFT )
		delete[] this->stripOrigins_NFFT;

	this->stripOrigins_NFFT = new UINTd[number_of_strips_NFFT];
	memcpy( stripOrigins_NFFT, pre_NFFT->stripOrigins_NFFT, number_of_strips_NFFT*sizeof(UINTd) );

	if( this->stripLengths_NFFT )
		delete[] this->stripLengths_NFFT;

	this->stripLengths_NFFT = new unsigned int[number_of_strips_NFFT];
	memcpy( stripLengths_NFFT, pre_NFFT->stripLengths_NFFT, number_of_strips_NFFT*sizeof(unsigned int) );

	unsigned int number_of_domains_NFFT_H = prod(pre_NFFT_H->domain_count_grid);

	if( this->domainsMap_NFFT_H )
		delete[] this->domainsMap_NFFT_H;

	this->domainsMap_NFFT_H = new unsigned int[number_of_domains_NFFT_H];
	memcpy( domainsMap_NFFT_H, pre_NFFT_H->domainsMap_NFFT_H, number_of_domains_NFFT_H*sizeof(unsigned int) );

//	this->number_of_threads_NFFT_H = pre_NFFT_H->number_of_threads_NFFT_H;
	this->number_of_strips_NFFT_H = pre_NFFT_H->number_of_strips_NFFT_H;

	if( this->stripsMap_NFFT_H )
		delete[] this->stripsMap_NFFT_H;

	this->stripsMap_NFFT_H = new uint2[number_of_domains_NFFT_H];
	memcpy( stripsMap_NFFT_H, pre_NFFT_H->stripsMap_NFFT_H, number_of_domains_NFFT_H*sizeof(uint2) );

	if( this->strips_NFFT_H )
		delete[] this->strips_NFFT_H;

	this->strips_NFFT_H = new uint2[number_of_strips_NFFT_H];
	memcpy( strips_NFFT_H, pre_NFFT_H->strips_NFFT_H, number_of_strips_NFFT_H*sizeof(uint2) );

	if( !sample_positions || !stripsMap_NFFT || !stripsDir_NFFT || !stripOrigins_NFFT || !stripLengths_NFFT || !domainsMap_NFFT_H || !stripsMap_NFFT_H || !strips_NFFT_H ){
		printf("\nOut of memory error! Quitting\n");
		exit(1);
	}

	this->max_strips_per_domain_NFFT = pre_NFFT->max_strips_per_domain_NFFT;

	this->domainsMapDevPtr_NFFT_H = pre_NFFT_H->domainsMapDevPtr_NFFT_H;
	this->stripsMapDevPtr_NFFT_H = pre_NFFT_H->stripsMapDevPtr_NFFT_H;
	this->stripsDevPtr_NFFT_H = pre_NFFT_H->stripsDevPtr_NFFT_H;

	this->stripsDirDevPtr_NFFT = pre_NFFT->stripsDirDevPtr_NFFT;
	this->stripOriginsDevPtr_NFFT = pre_NFFT->stripOriginsDevPtr_NFFT;
	this->stripsMapDevPtr_NFFT = pre_NFFT->stripsMapDevPtr_NFFT;
	this->stripLengthsDevPtr_NFFT = pre_NFFT->stripLengthsDevPtr_NFFT;

	delete pre_NFFT;
	delete pre_NFFT_H;

	this->successfully_preprocessed = true;

	return true;
}


	
template< class UINTd, class FLOATd, char TYPE >
NFFT_plan< UINTd, FLOATd, TYPE >::~NFFT_plan()
{
// Use (and implement) cleanup...
/*
	if( sample_positions ) delete[] sample_positions;
	if( domainsMap_NFFT ) delete[] domainsMap_NFFT;
	if( stripsMap_NFFT ) delete[] stripsMap_NFFT;
	if( stripsDir_NFFT ) delete[] stripsDir_NFFT;
	if( stripOrigins_NFFT ) delete[] stripOrigins_NFFT;
	if( stripLengths_NFFT ) delete[] stripLengths_NFFT;
*/
}

template< class UINTd, class FLOATd, char TYPE >
NFFT_H_plan< UINTd, FLOATd, TYPE >::~NFFT_H_plan()
{
// Use (and implement) cleanup...
/*
	if( sample_positions ) delete[] sample_positions;
	if( domainsMap_NFFT_H ) delete[] domainsMap_NFFT_H;
	if( stripsMap_NFFT_H ) delete[] stripsMap_NFFT_H;
	if( strips_NFFT_H ) delete[] strips_NFFT_H;
*/
}

template< class UINTd, class FLOATd, char TYPE >
NFFT_iteration_plan< UINTd, FLOATd, TYPE >::~NFFT_iteration_plan()
{
// Use (and implement) cleanup...
/*
	delete[] sample_positions;
	delete[] domainsMap_NFFT;
	delete[] stripsMap_NFFT;
	delete[] stripsDir_NFFT;
	delete[] stripOrigins_NFFT;
	delete[] stripLengths_NFFT;
	delete[] domainsMap_NFFT_H;
	delete[] stripsMap_NFFT_H;
	delete[] strips_NFFT_H;
*/
}

sample_strip::sample_strip( unsigned int idx )
{
	startIndex = idx;
	count = 1;
}


gridDomain::gridDomain()
{ 
#ifdef CIRCULAR_KERNEL
	origin = 0x0;
	squared_radius = 0.0f;
#endif
}

gridDomain::gridDomain( RealFloatArray &dimensions )
{ 
#ifdef CIRCULAR_KERNEL
	origin = 0x0;
	squared_radius = 0.0f;
#endif
	set( dimensions );
}

void 
gridDomain::set( RealFloatArray &dimensions )
{ 
#ifdef _DEBUG
	if (dimensions.get_number_of_dimensions() != 2)
	{
		cout << "Invalid array format in <gridDomain>. Number of dimensions MUST be 2." << endl;
		return;
	}

	if (dimensions.get_size(0) != 2)
	{
		cout << "Invalid array format in <gridDomain>. Width MUST be 2." << endl;
		return;
	}
#endif

	this->dimensions = dimensions;
}
	
#ifdef CIRCULAR_KERNEL
void 
gridDomain::set( unsigned int *origin, float squared_radius )
{ 
	this->squared_radius = squared_radius;
	this->origin = (unsigned int*) malloc( dimensions.get_size(1)*sizeof(unsigned int) );
	for( int i=0; i<dimensions.get_size(1); i++ )
		this->origin[i] = origin[i];
}
#endif

RealFloatArray&
gridDomain::get()
{ 
	return dimensions;
}

bool 
gridDomain::inside( double *pt )
{
	for( int i=0; i<dimensions.get_size(1); i++ )
		if( pt[i]<dimensions.get_data_ptr()[2*i+0] || pt[i]>dimensions.get_data_ptr()[2*i+1] )
			return false;

	return true;
}

bool 
gridDomain::inside( float *pt )
{
	for( int i=0; i<dimensions.get_size(1); i++ )
		if( pt[i]<dimensions.get_data_ptr()[2*i+0] || pt[i]>dimensions.get_data_ptr()[2*i+1] )
			return false;

	return true;
}

bool 
gridDomain::inside( unsigned int *pt )
{
	for( int i=0; i<dimensions.get_size(1); i++ )
		if( pt[i]<dimensions.get_data_ptr()[2*i+0] || pt[i]>dimensions.get_data_ptr()[2*i+1] )
			return false;

	return true;
}

#ifdef CIRCULAR_KERNEL
bool 
gridDomain::inside_circle( float *pt )
{
	float squared_dist = 0.0f;
	for( int i=0; i<dimensions.get_size(1); i++ )
		squared_dist += ((pt[i]-(float)origin[i])*(pt[i]-(float)origin[i]));

	return (squared_dist<=squared_radius);
}
#endif

extern "C" NFFT_H_plan< uint2, float2, 0>* 
preprocess_deapodization_2D( uint2 matrix_size, uint2 matrix_size_os, float W, uint2 domain_size_grid )
{
	// Generate trajectory with a single point at the origin
	RealFloatArray trajectory( 1,2 );
	trajectory.get_data_ptr()[0] = 0.0f; trajectory.get_data_ptr()[1] = 0.0f;

	uint2 fixed_dims = make_uint2( 0,0 );
	
	NFFT_H_plan< uint2, float2, 0> *plan = (NFFT_H_plan< uint2, float2, 0>*) new NFFT_H_plan< uint2, float2, 0>( *((uint2*)&matrix_size), *((uint2*)&matrix_size_os), fixed_dims, *((uint2*)&domain_size_grid), 1, W  );
	plan->preprocess( &trajectory );

	return plan;
}

extern "C" NFFT_H_plan< uint3, float3, 0>* 
preprocess_deapodization_3D( uint3 matrix_size, uint3 matrix_size_os, float W, uint3 domain_size_grid )
{
	// Generate trajectory with a single point at the origin
	RealFloatArray trajectory( 1,3 );
	trajectory.get_data_ptr()[0] = 0.0f; trajectory.get_data_ptr()[1] = 0.0f; trajectory.get_data_ptr()[2] = 0.0f;

	uint3 fixed_dims = make_uint3( 0,0,0 );

	NFFT_H_plan< uint3, float3, 0> *plan = (NFFT_H_plan< uint3, float3, 0>*) new NFFT_H_plan< uint3, float3, 0>( *((uint3*)&matrix_size), *((uint3*)&matrix_size_os), fixed_dims, *((uint3*)&domain_size_grid), 1, W );
	plan->preprocess( &trajectory );

	return plan;
}

extern "C" NFFT_H_plan< uint4, float4, 0>* 
preprocess_deapodization_4D( uint4 matrix_size, uint4 matrix_size_os, float W, uint4 domain_size_grid )
{
	// Generate trajectory with a single point at the origin
	RealFloatArray trajectory( 1,2 );
	trajectory.get_data_ptr()[0] = 0.0f; trajectory.get_data_ptr()[1] = 0.0f;
	trajectory.get_data_ptr()[2] = 0.0f; trajectory.get_data_ptr()[3] = 0.0f;

	uint4 fixed_dims = make_uint4( 0,0,0,0 );

	NFFT_H_plan< uint4, float4, 0> *plan = (NFFT_H_plan< uint4, float4, 0>*) new NFFT_H_plan< uint4, float4, 0>( *((uint4*)&matrix_size), *((uint4*)&matrix_size_os), fixed_dims, *((uint4*)&domain_size_grid), 1, W );
	plan->preprocess( &trajectory );

	return plan;
}

extern "C" void 
preprocess_deapodization_2D_cleanup( mr_recon::NFFT_H_plan< uint2, float2, 0> *pre )
{
	delete pre;
}

extern "C" void
preprocess_deapodization_3D_cleanup( mr_recon::NFFT_H_plan< uint3, float3, 0> *pre )
{
	delete pre;
}

extern "C" void
preprocess_deapodization_4D_cleanup( mr_recon::NFFT_H_plan< uint4, float4, 0> *pre )
{
	delete pre;
}


// Template instantiation

template class NFFT_plan< uint2, float2, 0>;
template class NFFT_plan< uint3, float3, 0>;
template class NFFT_plan< uint4, float4, 0>;

template class NFFT_H_plan< uint2, float2, 0>;
template class NFFT_H_plan< uint3, float3, 0>;
template class NFFT_H_plan< uint4, float4, 0>;

template class NFFT_iteration_plan< uint2, float2, 0>;
template class NFFT_iteration_plan< uint3, float3, 0>;
template class NFFT_iteration_plan< uint4, float4, 0>;

template class NFFT_plan< uint2, float2, 1>;
template class NFFT_plan< uint3, float3, 1>;
template class NFFT_plan< uint4, float4, 1>;

template class NFFT_H_plan< uint2, float2, 1>;
template class NFFT_H_plan< uint3, float3, 1>;
template class NFFT_H_plan< uint4, float4, 1>;

template class NFFT_iteration_plan< uint2, float2, 1>;
template class NFFT_iteration_plan< uint3, float3, 1>;
template class NFFT_iteration_plan< uint4, float4, 1>;

template class NFFT_plan< uint2, float2, 2>;
template class NFFT_plan< uint3, float3, 2>;
template class NFFT_plan< uint4, float4, 2>;

template class NFFT_H_plan< uint2, float2, 2>;
template class NFFT_H_plan< uint3, float3, 2>;
template class NFFT_H_plan< uint4, float4, 2>;

template class NFFT_iteration_plan< uint2, float2, 2>;
template class NFFT_iteration_plan< uint3, float3, 2>;
template class NFFT_iteration_plan< uint4, float4, 2>;

template mr_recon::NFFT_plan< uint2, float2, 0>* mr_recon::preprocess_NFFT< uint2, float2, 0>( uint2, uint2, uint2, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_plan< uint3, float3, 0>* mr_recon::preprocess_NFFT< uint3, float3, 0>( uint3, uint3, uint3, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_plan< uint4, float4, 0>* mr_recon::preprocess_NFFT< uint4, float4, 0>( uint4, uint4, uint4, unsigned int, unsigned int, float, RealFloatArray * );

template mr_recon::NFFT_H_plan< uint2, float2, 0>* mr_recon::preprocess_NFFT< uint2, float2, 0>( uint2, uint2, uint2, uint2, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_H_plan< uint3, float3, 0>* mr_recon::preprocess_NFFT< uint3, float3, 0>( uint3, uint3, uint3, uint3, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_H_plan< uint4, float4, 0>* mr_recon::preprocess_NFFT< uint4, float4, 0>( uint4, uint4, uint4, uint4, unsigned int, float, RealFloatArray* );

template mr_recon::NFFT_iteration_plan< uint2, float2, 0>* mr_recon::preprocess_NFFT< uint2, float2, 0>( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint3, float3, 0>* mr_recon::preprocess_NFFT< uint3, float3, 0>( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint4, float4, 0>* mr_recon::preprocess_NFFT< uint4, float4, 0>( uint4, uint4, uint4, uint4, unsigned int, unsigned int, float, RealFloatArray* );

template mr_recon::NFFT_plan< uint2, float2, 1>* mr_recon::preprocess_NFFT< uint2, float2, 1>( uint2, uint2, uint2, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_plan< uint3, float3, 1>* mr_recon::preprocess_NFFT< uint3, float3, 1>( uint3, uint3, uint3, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_plan< uint4, float4, 1>* mr_recon::preprocess_NFFT< uint4, float4, 1>( uint4, uint4, uint4, unsigned int, unsigned int, float, RealFloatArray * );

template mr_recon::NFFT_H_plan< uint2, float2, 1>* mr_recon::preprocess_NFFT< uint2, float2, 1>( uint2, uint2, uint2, uint2, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_H_plan< uint3, float3, 1>* mr_recon::preprocess_NFFT< uint3, float3, 1>( uint3, uint3, uint3, uint3, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_H_plan< uint4, float4, 1>* mr_recon::preprocess_NFFT< uint4, float4, 1>( uint4, uint4, uint4, uint4, unsigned int, float, RealFloatArray* );

template mr_recon::NFFT_iteration_plan< uint2, float2, 1>* mr_recon::preprocess_NFFT< uint2, float2, 1>( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint3, float3, 1>* mr_recon::preprocess_NFFT< uint3, float3, 1>( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint4, float4, 1>* mr_recon::preprocess_NFFT< uint4, float4, 1>( uint4, uint4, uint4, uint4, unsigned int, unsigned int, float, RealFloatArray* );

template mr_recon::NFFT_plan< uint2, float2, 2>* mr_recon::preprocess_NFFT< uint2, float2, 2>( uint2, uint2, uint2, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_plan< uint3, float3, 2>* mr_recon::preprocess_NFFT< uint3, float3, 2>( uint3, uint3, uint3, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_plan< uint4, float4, 2>* mr_recon::preprocess_NFFT< uint4, float4, 2>( uint4, uint4, uint4, unsigned int, unsigned int, float, RealFloatArray * );

template mr_recon::NFFT_H_plan< uint2, float2, 2>* mr_recon::preprocess_NFFT< uint2, float2, 2>( uint2, uint2, uint2, uint2, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_H_plan< uint3, float3, 2>* mr_recon::preprocess_NFFT< uint3, float3, 2>( uint3, uint3, uint3, uint3, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_H_plan< uint4, float4, 2>* mr_recon::preprocess_NFFT< uint4, float4, 2>( uint4, uint4, uint4, uint4, unsigned int, float, RealFloatArray* );

template mr_recon::NFFT_iteration_plan< uint2, float2, 2>* mr_recon::preprocess_NFFT< uint2, float2, 2>( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint3, float3, 2>* mr_recon::preprocess_NFFT< uint3, float3, 2>( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint4, float4, 2>* mr_recon::preprocess_NFFT< uint4, float4, 2>( uint4, uint4, uint4, uint4, unsigned int, unsigned int, float, RealFloatArray* );

