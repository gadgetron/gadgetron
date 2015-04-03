/*
 * osSenseOperator.h
 *
 *  Created on: Mar 31, 2015
 *      Author: dch
 */

#pragma once
#include "subsetOperator.h"
#include <boost/shared_ptr.hpp>

namespace Gadgetron {

template<class ARRAY, class SENSE > class osSenseOperator: public Gadgetron::subsetOperator<ARRAY>, public FFTOperator {
public:
	osSenseOperator(){};
	virtual ~osSenseOperator(){};

	virtual void set_csm(boos::shared_ptr<ARRAY> csm){
		this->csm = csm;
	}
	virtual boost::shared_ptr<ARRAY> get_csm(){
		return csm;
	}
	virtual int get_number_of_subsets() override {
		return csm->size(csm->get_number_of_dimensions()-1);
	}

	virtual void mult_M(ARRAY* in, ARRAY* out, int subset, bool accumulate) override {
		auto subsize = *csm->get_dimensions();
		subsize.pop_back();
		auto num_elements = std::accumulate(subsize.begin(),subsize.end(),1,std::multiplies<typename ARRAY::element_type>());
		ARRAY sub_csm(subsize,csm->get_data_ptr()+num_elements*subset);
		ARRAY tmp(in->get_dimensions());
		csm_mult_M(in,&tmp,&sub_csm);
		FFTOperator::mult_M(&tmp,out,accumulate);
	}

	virtual void mult_MH(ARRAY* in, ARRAY* out, int subset, bool accumulate) override {
		auto subsize = *csm->get_dimensions();
		subsize.pop_back();
		auto num_elements = std::accumulate(subsize.begin(),subsize.end(),1,std::multiplies<typename ARRAY::element_type>());
		ARRAY sub_csm(subsize,csm->get_data_ptr()+num_elements*subset);

		ARRAY tmp(out->get_dimensions());

		FFTOperator::mult_MH(in,&tmp,accumulate);
		csm_mult_M(&tmp,out,&sub_csm);

	}



protected:
	boost::shared_ptr<ARRAY> csm;
};

} /* namespace Gadgetron */

