/** \file multiplicationOperatorContainer.h
    \brief Operator used to chain together (concatenate) a series of operators by multiplication.
 */

#pragma once

#include "linearOperator.h"
#include <iostream>
#include <vector>

namespace Gadgetron {

template <class ARRAY_TYPE> class multiplicationOperatorContainer : public linearOperator<ARRAY_TYPE> {
    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;

  public:
    multiplicationOperatorContainer() : linearOperator<ARRAY_TYPE>() {}
    virtual ~multiplicationOperatorContainer() {}

    // Set/get domain and codomain dimensions.
    //

    virtual void set_domain_dimensions(const std::vector<unsigned int>&) {
        throw std::runtime_error("Warning: multiplicationOperatorContainer::set_domain_dimensions : dimensions "
                                 "ignored, using dimensions of the individual operators instead");
    }

    virtual void set_codomain_dimensions(const std::vector<size_t>&) {
        throw std::runtime_error("Warning: multiplicationOperatorContainer::set_codomain_dimensions : dimensions "
                                 "ignored, using dimensions of the individual operators instead");
    }

    virtual std::vector<size_t> get_domain_dimensions() {
        if (operators_.size() == 0)
            return std::vector<size_t>();
        else
            return operators_[0]->get_domain_dimensions();
    }

    virtual std::vector<size_t> get_codomain_dimensions() {
        if (operators_.size() == 0)
            return std::vector<size_t>();
        else
            return operators_[operators_.size() - 1]->get_codomain_dimensions();
    }

    virtual void set_weight(REAL weight) {
        REAL op_weight = REAL(1);
        for (int i = 0; i < operators_.size(); i++)
            op_weight *= operators_[i]->get_weight();
        this->weight_ = weight * op_weight;
    }

    // Add operator to the container
    //
    void add_operator(boost::shared_ptr<linearOperator<ARRAY_TYPE>> op) {
        if (op.get() == 0x0) {
            throw std::runtime_error("Error: multiplicationOperatorContainer::add_operator : illegal operator");
        }

        // All operators needs the domain and codomain dimensions set
        //
        if (op->get_domain_dimensions().empty()) {
            throw std::runtime_error(
                "Error: multiplicationOperatorContainer::add_operator : domain dimensions not set on operator");
        }
        if (op->get_codomain_dimensions().empty()) {
            throw std::runtime_error(
                "Error: multiplicationOperatorContainer::add_operator : codomain dimensions not set on operator");
        }
        operators_.push_back(op);
        this->weight_ *= op->get_weight();
    }

    virtual void mult_M(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) {
        if (operators_.size() == 0) {
            throw std::runtime_error("Error: multiplicationOperatorContainer::mult_M : no operators added");
        }

        ARRAY_TYPE *tmp_in = in, *tmp_out = 0x0;
        ARRAY_TYPE ping, pong;

        if (operators_.size() > 1) {
            ping.create(operators_[0]->get_codomain_dimensions());
            tmp_out = &ping;
        } else {
            tmp_out = out;
        }

        // Loop over operators
        //
        for (int i = 0; i < operators_.size(); i++) {

            operators_[i]->mult_M(tmp_in, tmp_out, (i == operators_.size() - 1) ? accumulate : false);

            ARRAY_TYPE* tmp_tmp_out = (i == 0) ? &pong : tmp_in;
            tmp_in = tmp_out;

            if (operators_.size() > 2 && i < operators_.size() - 2) {
                tmp_tmp_out->create(operators_[i + 1]->get_codomain_dimensions());
                tmp_out = tmp_tmp_out;
            } else if (i == operators_.size() - 2) {
                tmp_out = out;
            }
        }
    }

    virtual void mult_MH(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) {
        if (operators_.size() == 0) {
            throw std::runtime_error("Error: multiplicationOperatorContainer::mult_MH : no operators added");
        }

        ARRAY_TYPE *tmp_in = in, *tmp_out = 0x0;
        ARRAY_TYPE ping, pong;

        if (operators_.size() > 1) {
            ping.create(operators_[operators_.size() - 1]->get_domain_dimensions());
            tmp_out = &ping;
        } else {
            tmp_out = out;
        }

        // Loop over operators
        //
        for (int i = operators_.size() - 1; i >= 0; i--) {

            operators_[i]->mult_MH(tmp_in, tmp_out, (i == 0) ? accumulate : false);

            ARRAY_TYPE* tmp_tmp_out = (i == operators_.size() - 1) ? &pong : tmp_in;
            tmp_in = tmp_out;

            if (i > 1) {
                tmp_tmp_out->create(operators_[i - 1]->get_domain_dimensions());
                tmp_out = tmp_tmp_out;
            } else if (i == 1) {
                tmp_out = out;
            }
        }
    }

    virtual boost::shared_ptr<linearOperator<ARRAY_TYPE>> clone() { return linearOperator<ARRAY_TYPE>::clone(this); }

  protected:
    std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE>>> operators_;
};
} // namespace Gadgetron
