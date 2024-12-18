#pragma once

#include "cuNDArray.h"
#include "complext.h"

namespace Gadgetron{

// Multiply with coil sensitivities
//

template< class REAL, unsigned int D> void
csm_mult_M( cuNDArray< complext<REAL> > *in,
	    cuNDArray< complext<REAL> > *out,
	    cuNDArray< complext<REAL> > *csm );


// Multiply with adjoint of coil sensitivities
//

template< class REAL, unsigned int D> void
csm_mult_MH( cuNDArray< complext<REAL> > *in,
	     cuNDArray< complext<REAL> > *out,
	     cuNDArray< complext<REAL> > *csm );
}
