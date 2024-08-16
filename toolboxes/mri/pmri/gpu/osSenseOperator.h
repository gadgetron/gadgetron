/*
 * osSenseOperator.h
 *
 *  Created on: Mar 31, 2015
 *      Author: dch
 */

#pragma once
#include "subsetOperator.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "sense_utilities.h"
namespace Gadgetron {

template<class ARRAY, unsigned int D, class FFTOperator> class osSenseOperator: public virtual Gadgetron::subsetOperator<ARRAY>, public virtual FFTOperator {
public:
	osSenseOperator(): coils_per_subset(1){};
	virtual ~osSenseOperator(){};
	typedef typename ARRAY::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::Type REAL;

	virtual void set_csm(boost::shared_ptr<ARRAY> csm){
		if (csm->get_size(csm->get_number_of_dimensions()-1)%coils_per_subset != 0)
			throw std::runtime_error("osSenseOperator: number of coils in coil sensitivity map must be divisible by coils per subset");
		this->csm = csm;
	}
	virtual boost::shared_ptr<ARRAY> get_csm(){
		return csm;
	}
	virtual int get_number_of_subsets() override {
		return csm->get_size(csm->get_number_of_dimensions()-1)/coils_per_subset;
	}

	virtual void set_coils_per_subset(unsigned int coils_per_subset){
		this->coils_per_subset = coils_per_subset;
	}

	virtual void mult_M(ARRAY* in, ARRAY* out, int subset, bool accumulate) override {
		auto subsize = csm->get_dimensions();
		subsize.pop_back();
		subsize.push_back(coils_per_subset);
		auto num_elements = std::accumulate(subsize.begin(),subsize.end(),1,std::multiplies<size_t>());
		ARRAY sub_csm(subsize,csm->get_data_ptr()+num_elements*subset);
		auto in_dims = in->get_dimensions();
		in_dims.push_back(coils_per_subset);
		ARRAY tmp(in_dims);
		csm_mult_M<REAL,D>(in,&tmp,&sub_csm);
		FFTOperator::mult_M(&tmp,out,accumulate);
	}

	virtual void mult_M(ARRAY* in, ARRAY* out, bool accumulate) final {
		subsetOperator<ARRAY>::mult_M(in,out,accumulate);
	}

	virtual void mult_MH(ARRAY* in, ARRAY* out, bool accumulate) final {
		subsetOperator<ARRAY>::mult_MH(in,out,accumulate);
	}

	virtual void mult_MH_M(ARRAY* in, ARRAY* out, bool accumulate) final {

		ARRAY tmp(this->codomain_dims_);
		mult_M(in,&tmp,false);
		mult_MH(&tmp,out,accumulate);
	}

	virtual void mult_MH(ARRAY* in, ARRAY* out, int subset, bool accumulate) override {
		auto subsize = csm->get_dimensions();
		subsize.pop_back();
		subsize.push_back(coils_per_subset);
		auto num_elements = std::accumulate(subsize.begin(),subsize.end(),1,std::multiplies<size_t>());
		ARRAY sub_csm(subsize,csm->get_data_ptr()+num_elements*subset);

		auto out_dims = out->get_dimensions();
		out_dims.push_back(coils_per_subset);
		ARRAY tmp(out_dims);

		FFTOperator::mult_MH(in,&tmp,accumulate);

		csm_mult_MH<REAL,D>(&tmp,out,&sub_csm);

	}

	virtual std::vector<size_t> get_codomain_dimensions(int subset) override{
		auto codom_dims = this->codomain_dims_;
		codom_dims.back() = coils_per_subset;
		return codom_dims;
	}



protected:
	boost::shared_ptr<ARRAY> csm;
	unsigned int coils_per_subset;
};

} /* namespace Gadgetron */

