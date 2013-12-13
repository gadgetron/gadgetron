#pragma once

#include "linearOperator.h"
#include <numeric>
#include <functional>
namespace Gadgetron{


template<class ARRAY_TYPE> class subsetOperator : public linearOperator<ARRAY_TYPE>{
private:
  typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
  typedef typename realType<ELEMENT_TYPE>::Type REAL;
public:

	subsetOperator(int _number_of_subsets) : number_of_subsets(_number_of_subsets),
	subset_dimensions(number_of_subsets,std::vector<size_t>()){};

	virtual ~subsetOperator(){};
	virtual void mult_M(ARRAY_TYPE* in, ARRAY_TYPE* out, int subset, bool accumulate)=0;
	virtual void mult_MH(ARRAY_TYPE* in, ARRAY_TYPE* out, int subset, bool accumulate)=0;
	virtual void mult_MH_M(ARRAY_TYPE* in, ARRAY_TYPE* out, int subset, bool accumulate)=0;

	virtual void mult_M(ARRAY_TYPE* in, ARRAY_TYPE* out,bool accumulate){
		if (!accumulate) clear(out);
		std::vector<boost::shared_ptr<ARRAY_TYPE> > projections = projection_subsets(out);

		for (int i = 0; i < number_of_subsets; i++) mult_M(in,projections[i].get(),i,true);
	}

	virtual void mult_MH(ARRAY_TYPE* in, ARRAY_TYPE* out,bool accumulate){
			if (!accumulate) clear(out);
			std::vector<boost::shared_ptr<ARRAY_TYPE> > projections = projection_subsets(in);
			for (int i = 0; i < number_of_subsets; i++) mult_MH(projections[i].get(),out,i,true);
	}

	virtual boost::shared_ptr< std::vector<size_t> > get_codomain_dimensions(int subset){
		return boost::shared_ptr< std::vector<size_t> >(new std::vector<size_t>(subset_dimensions[subset]));
	}
/*
 	virtual void set_codomain_subsets(std::vector< std::vector<unsigned int> > & _dims){
		codomain_dimensions = std::vector< std::vector<unsigned int> >(_dims);
	}
*/
	int get_number_of_subsets(){return number_of_subsets;}


	virtual std::vector<boost::shared_ptr<ARRAY_TYPE> > projection_subsets(ARRAY_TYPE* projections){

		std::vector<boost::shared_ptr<ARRAY_TYPE> > res;
		ELEMENT_TYPE* curPtr = projections->get_data_ptr();
		for (int subset = 0; subset < number_of_subsets; subset++){
			res.push_back(boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE(&subset_dimensions[subset],curPtr)));
			curPtr += std::accumulate(subset_dimensions[subset].begin(),subset_dimensions[subset].end(),1,std::multiplies<unsigned int>());
		}
		return res;
	}

protected:
	int number_of_subsets;
	std::vector< std::vector<size_t> > subset_dimensions;
};
}
