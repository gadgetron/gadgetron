#include "cuNDArray.h"
#include <thrust/iterator/zip_iterator.h>
#include "vector_td.h"
#include "vector_td_utilities.h"
#include "trajectory_utils.h"
using namespace Gadgetron;

struct traj_filter_functor{

	traj_filter_functor(float limit){
		limit_ = limit;
	}
	__device__  thrust::tuple<floatd2,float> operator()(thrust::tuple<floatd2,float> tup){

		floatd2 traj = thrust::get<0>(tup);
		float dcw = thrust::get<1>(tup);
		if ( abs(traj[0]) > limit_ || abs(traj[1]) > limit_)
			dcw = 0;

		return thrust::tuple<floatd2,float>(traj,dcw);

	}

	float limit_;
};

boost::shared_ptr< cuNDArray<float> > Gadgetron::filter_dcw(cuNDArray<floatd2>* traj, cuNDArray<float>* dcw,float limit){
	cuNDArray<float>* dcw_new = new cuNDArray<float>(*dcw);

	thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(traj->begin(),dcw_new->begin())),
			thrust::make_zip_iterator(thrust::make_tuple(traj->end(),dcw_new->end())),
			thrust::make_zip_iterator(thrust::make_tuple(traj->begin(),dcw_new->begin())),
			traj_filter_functor(limit));

	return boost::shared_ptr<cuNDArray<float> >(dcw_new);


}
