#include "linearOperator.h"
#include "hoCuNDArray.h"
#include <vector>
#include "GPUTimer.h"
#pragma once

namespace Gadgetron{
template<class REAL>
class hoCuOperatorPathBackprojection : public linearOperator<hoCuNDArray<REAL> > {
 public:
 hoCuOperatorPathBackprojection() : linearOperator<hoCuNDArray<REAL> >() {
	 rescale_dirs=true;
    }
    virtual ~hoCuOperatorPathBackprojection() {}

    virtual void mult_M( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate = false );
    virtual void mult_MH( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate = false );
    virtual void mult_MH_M( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate = false );

    virtual void setup(boost::shared_ptr< hoCuNDArray<vector_td<REAL,3> > > splines, vector_td<REAL,3> physical_dims,boost::shared_ptr< hoCuNDArray< REAL > > projections, vector_td<REAL,3> origin, REAL background=0);
    virtual void setup(boost::shared_ptr< hoCuNDArray<vector_td<REAL,3> > > splines, vector_td<REAL,3> physical_dims, boost::shared_ptr< hoCuNDArray< REAL > > projections, boost::shared_ptr< hoCuNDArray<REAL> > weights, vector_td<REAL,3> origin, REAL background = 0);

    virtual boost::shared_ptr< linearOperator<hoCuNDArray<REAL> > > clone() {
       return linearOperator<hoCuNDArray<REAL> >::clone(this);
     }
    void rescale_directions(bool val){rescale_dirs=val;}
 protected:
    boost::shared_ptr< hoCuNDArray<vector_td<REAL,3> > > splines;

    boost::shared_ptr< hoCuNDArray<REAL> > weights;
    vector_td<REAL,3> physical_dims;
    vector_td<REAL,3> origin;

    REAL background;

    bool rescale_dirs;
    size_t calculate_batch_size();

};
}
