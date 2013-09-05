#pragma once

namespace Gadgetron{

  template<class ARRAY_TYPE, class TV_OPERATOR, class REAL> class tvPicsOperator 
    : public generalOperator<ARRAY_TYPE>
  {
  public:
    
    tvPicsOperator() : generalOperator<ARRAY_TYPE>() {}
    virtual ~tvPicsOperator() {}

    void set_prior(boost::shared_ptr<ARRAY_TYPE> prior){
      prior_ = prior;
    }

    virtual void gradient(ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate){
      ARRAY_TYPE tmp = *in;
      tmp -= *prior_;
      op_.gradient(&tmp, out, accumulate);
    }

    virtual REAL magnitude(ARRAY_TYPE *x){
    	ARRAY_TYPE tmp = *x;
    	tmp -= *prior_;
    	return op_.magnitude(&tmp);
    }
    void set_limit(REAL limit){
      op_.set_limit(limit);
    }

    virtual void set_weight(REAL weight){
      op_.set_weight(weight);
    }

    virtual REAL get_weight(){
      return op_.get_weight();
    }

  protected:
    TV_OPERATOR op_;
    boost::shared_ptr<ARRAY_TYPE> prior_;
  };
}
