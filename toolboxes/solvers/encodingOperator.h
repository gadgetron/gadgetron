#pragma once

#include "linearOperator.h"
#include <iostream>
#include <vector>
#include <boost/smart_ptr.hpp>


template <class REAL, class ARRAY_TYPE> class encodingOperator
	: public linearOperator<REAL, ARRAY_TYPE>
{
	protected:
		std::vector< boost::shared_ptr<linearOperator<REAL, ARRAY_TYPE> > > operators_;
		std::vector<unsigned int> offsets_;

	public:
		encodingOperator() : linearOperator<REAL,ARRAY_TYPE>() {}
		virtual ~encodingOperator(){}
		inline std::vector<unsigned int> get_domain_dimensions(int i) { return operators_[i]->get_domain_dimensions(); }
		inline std::vector<unsigned int> get_codomain_dimensions(int i) { return operators_[i]->get_codomain_dimensions(); }
		boost::shared_ptr< ARRAY_TYPE> create_codomain() {
			ARRAY_TYPE* codomain = new ARRAY_TYPE;
			codomain->create(this->get_codomain_dimensions());
			return boost::shared_ptr<ARRAY_TYPE>(codomain);
		}
		unsigned int get_offset(int i){
			return offsets_[i];
		}
		boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE > > get_operator(int i){
			return operators_[i];
		}

		void add_operator(boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE > > op){
			int elements = 1;
			std::vector<int> codomain= operators_.back()->get_codomain_dimensions();
			for (int i =0; i < codomain.size(); i++){
				elements *= codomain[i];
			}

			if (offsets_.size() == 0){
				offsets_.push_back(0);
			} else{


				offsets_.push_back(offsets_.back()+elements);
			}
			operators_.push_back(op);
			codomain_dims_[0]+= elements;
		}

		  virtual int mult_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false){
			  for (int i =0; i < operators_.size(); i++){
				  ARRAY_TYPE tmp;
				  tmp.create(operators_[i]->get_codomain_dimensions(),out->get_data_ptr()+offsets_[i]);
				  int res = operators_[i]->mult_M(in,&tmp,accumulate);
				  if (res < 0) return res;

			  }
			  return 0;
		  }

		  virtual int mult_MH( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false){
			  for (int i =0; i < operators_.size(); i++){
				  ARRAY_TYPE tmp;
				   tmp.create(operators_[i]->get_codomain_dimensions(),in->get_data_ptr()+offsets_[i]);
				   int res;
				   if (!accumulate && i == 0){
					   res=operators_[i]->mult_MH(&tmp,out,false);
				   } else {
					   res=operators_[i]->mult_MH(&tmp,out,true);
				   }
				   if (res < 0) return res;

			  }
			  return 0;
		  }
		  virtual int mult_MH_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false){
			  for (int i =0; i < operators_.size(); i++){
				  int res;
				   if (!accumulate && i == 0){
					   res = operators_[i]->mult_MH(in,out,false);
				   } else {
					   res= operators_[i]->mult_MH(in,out,true);
				   }
				   if (res < 0) return res;

			  }
			  return 0;
		  }


};
