/** \file   hoImageRegParametricTransformation.h
    \brief  Define the base class for the parametric geometry transformation in gadgetron registration
    \author Hui Xue
*/

#ifndef hoImageRegParametricTransformation_H_
#define hoImageRegParametricTransformation_H_

#pragma once

#include "hoImageRegTransformation.h"

namespace Gadgetron {

    /// parametric transformation, e.g. rigid and affine transformation or Free-Form Deformation
    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    class hoImageRegParametricTransformation : public hoImageRegTransformation<ValueType, DIn, DOut>
    {
    public:

        typedef hoImageRegTransformation<ValueType, DIn, DOut> BaseClass;
        typedef hoImageRegParametricTransformation<ValueType, DIn, DOut> Self;

        typedef ValueType T;

        typedef typename BaseClass::input_point_type input_point_type;
        typedef typename BaseClass::output_point_type output_point_type;

        typedef typename BaseClass::jacobian_parameter_type jacobian_parameter_type;
        typedef typename BaseClass::jacobian_position_type jacobian_position_type;

        /// every parameter can be active or inactive
        /// if inactive, this parameter will not be changed during optimization
        typedef enum { Inactive=0, Active, Unknown } ParaStatus;
        typedef std::vector<ParaStatus> ParaStatusType;

        hoImageRegParametricTransformation() : num_parameters_(0), BaseClass() {}
        virtual ~hoImageRegParametricTransformation() {}

        size_t get_number_of_parameters() const { return num_parameters_; }
        void set_number_of_parameters(size_t num) { num_parameters_ = num; para_status_.resize(num, Active); }

        // get/set the ith parameter
        virtual ValueType get_parameter(size_t i) const = 0;
        virtual void set_parameter(size_t i, ValueType v) = 0;

        ParaStatus get_para_status(size_t i) { GADGET_CHECK_THROW(i<num_parameters_); return this->para_status_[i]; }
        void set_para_status(size_t i, ParaStatus status) { GADGET_CHECK_THROW(i<num_parameters_); para_status_[i] = status; }

        virtual bool invertTransformation() = 0;

        virtual bool setIdentity() = 0;

        virtual bool transform(const T* pt_in, T* pt_out) const = 0;

        virtual bool transform(const T& xi, const T& yi, T& xo, T& yo) const = 0;

        virtual bool transform(const T& xi, const T& yi, const T& zi, T& xo, T& yo, T& zo) const = 0;

        virtual bool transform(const size_t* pt_in, T* pt_out) const = 0;
        virtual bool transform(const size_t* pt_in, size_t N, T* pt_out) const = 0;
        virtual bool transform(const size_t& xi, const size_t& yi, T& xo, T& yo) const = 0;
        virtual bool transform(const size_t* xi, const size_t* yi, size_t N, T* xo, T* yo) const = 0;
        virtual bool transform(const size_t& xi, const size_t& yi, const size_t& zi, T& xo, T& yo, T& zo) const = 0;
        virtual bool transform(const size_t* xi, const size_t* yi, const size_t* zi, size_t N, T* xo, T* yo, T* zo) const = 0;

        /// adjust transformation for the resolution pyramid, if the image coordinate is used
        /// sourceI2W and targetI2W: source and target image to world transformation matrix
        virtual bool adjustForResolutionPyramid(const hoMatrix<ValueType>& sourceI2W, const hoMatrix<ValueType>& targetI2W)
        {
            /// by default, the transformation is not changed
            return true;
        }

        /// compute jacobian matrix to parameters
        /// DOut*num_parameters_ matrix
        virtual bool jacobianParameter(const input_point_type& /*pos*/, jacobian_parameter_type& jac)
        {
            jac.createMatrix(DOut, num_parameters_);
            jac.setIdentity();
            return true;
        }

        /// compute jacobian matrix to spatial position
        /// DOut*DIn matrix
        virtual bool jacobianPosition(const input_point_type& /*pos*/, jacobian_position_type& jac)
        {
            jac.createMatrix(DOut, DIn);
            jac.setIdentity();
            return true;
        }

        /// serialize/deserialize the transformation
        virtual bool serialize(char*& buf, size_t& len) const;
        virtual bool deserialize(char* buf, size_t& len);

        virtual void print(std::ostream& os) const
        {
            using namespace std;
            os << "--------------Gagdgetron parametric geometry transformation -------------" << endl;
            os << "Input dimension is : " << DIn << endl;
            os << "Output dimension is : " << DOut << endl;

            std::string elemTypeName = std::string(typeid(T).name());
            os << "Transformation data type is : " << elemTypeName << std::endl;
            os << "Number of parameters is : " << num_parameters_ << endl;

            size_t i;
            os << "Status of parameters: " << endl;
            for ( i=0; i<this->num_parameters_; i++ )
            {
                os << "Para " << i << " : \t";
                if ( para_status_[i] == Active )
                {
                    os << "Active";
                }
                else if ( para_status_[i] == Inactive )
                {
                    os << "Inactive";
                }
                else
                {
                    os << "Unknown";
                }
                os << endl;
            }
        }

        virtual void printTransform(std::ostream& os) const
        {
            using namespace std;

            size_t i;
            size_t maxNum = 12;

            if ( this->num_parameters_< maxNum )
            {
                os << "[ ";
                for ( i=0; i<this->num_parameters_; i++ )
                {
                    os << this->get_parameter(i) << " \t";
                }
                os << " ]" << endl;
            }
            else
            {
                os << "[ ";
                for ( i=0; i<maxNum; i++ )
                {
                    os << this->get_parameter(i) << " \t";
                }
                os << " ... ]" << endl;
            }
        }

        virtual std::string transformationName() const
        {
            return std::string("hoImageRegParametricTransformation"); 
        }

        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

    protected:

        size_t num_parameters_;

        ParaStatusType para_status_;
    };

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    bool hoImageRegParametricTransformation<ValueType, DIn, DOut>::serialize(char*& buf, size_t& len) const 
    {
        try
        {
            if ( buf != NULL ) delete[] buf;

            size_t numOfPara = this->get_number_of_parameters();
            size_t totalLen = sizeof(ValueType)*numOfPara;

            buf = new char[totalLen];
            GADGET_CHECK_RETURN_FALSE(buf!=NULL);

            ValueType currPara;
            size_t ii, offset(0);
            for ( ii=0; ii<numOfPara; ii++ )
            {
                currPara = this->get_parameter(ii);
                memcpy(buf+offset, &currPara, sizeof(ValueType));
                offset += sizeof(ValueType);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegParametricTransformation<ValueType, DIn, DOut>::serialize(char*& buf, size_t& len) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    bool hoImageRegParametricTransformation<ValueType, DIn, DOut>::deserialize(char* buf, size_t& len)
    {
        try
        {
            size_t numOfPara = this->get_number_of_parameters();

            ValueType currPara;
            size_t ii, offset(0);
            for ( ii=0; ii<numOfPara; ii++ )
            {
                memcpy(&currPara, buf+offset, sizeof(ValueType));
                offset += sizeof(ValueType);
                this->set_parameter(ii, currPara);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegParametricTransformation<ValueType, DIn, DOut>::deserialize(char* buf, size_t& len) ... ");
            return false;
        }

        return true;
    }
}
#endif // hoImageRegParametricTransformation_H_
