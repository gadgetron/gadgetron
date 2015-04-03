/*
 * Base class for handling operations on a subset of the data. Is used as the operator class for all
 * ordered subset solvers.
 */
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

	subsetOperator(int _number_of_subsets) : number_of_subsets(_number_of_subsets){};
	subsetOperator() : number_of_subsets(1){};

	virtual ~subsetOperator(){};
	virtual void mult_M(ARRAY_TYPE* in, ARRAY_TYPE* out, int subset, bool accumulate)=0;
	virtual void mult_MH(ARRAY_TYPE* in, ARRAY_TYPE* out, int subset, bool accumulate)=0;
	virtual void mult_MH_M(ARRAY_TYPE* in, ARRAY_TYPE* out, int subset, bool accumulate)=0;

	virtual void mult_M(ARRAY_TYPE* in, ARRAY_TYPE* out,bool accumulate){
		if (!accumulate) clear(out);
		std::vector<boost::shared_ptr<ARRAY_TYPE> > projections = projection_subsets(out);

		for (int i = 0; i < this->get_number_of_subsets(); i++) mult_M(in,projections[i].get(),i,true);
	}

	virtual void mult_MH(ARRAY_TYPE* in, ARRAY_TYPE* out,bool accumulate){
			if (!accumulate) clear(out);
			std::vector<boost::shared_ptr<ARRAY_TYPE> > projections = projection_subsets(in);
			for (int i = 0; i < this->get_number_of_subsets(); i++) mult_MH(projections[i].get(),out,i,true);
	}

	virtual boost::shared_ptr< std::vector<size_t> > get_codomain_dimensions(int subset)=0;
/*
 	virtual void set_codomain_subsets(std::vector< std::vector<unsigned int> > & _dims){
		codomain_dimensions = std::vector< std::vector<unsigned int> >(_dims);
	}
*/
	virtual int get_number_of_subsets(){return number_of_subsets;}



	virtual std::vector<boost::shared_ptr<ARRAY_TYPE> > projection_subsets(ARRAY_TYPE* projections){

		std::vector<boost::shared_ptr<ARRAY_TYPE> > res;
		ELEMENT_TYPE* curPtr = projections->get_data_ptr();
		for (int subset = 0; subset < this->get_number_of_subsets(); subset++){
			std::vector<size_t> subset_dim = *get_codomain_dimensions(subset);
			res.push_back(boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE(&subset_dim,curPtr)));
			curPtr += std::accumulate(subset_dim.begin(),subset_dim.end(),1,std::multiplies<unsigned int>());
		}
		return res;
	}
protected:
	virtual void set_number_of_subsets(int nsubsets){ number_of_subsets=nsubsets;}
private:
	int number_of_subsets;

};
}
