/** \file   hoImageRegHomogenousTransformation.h
    \brief  Define the class for the homogenous geometry transformation in gadgetron registration
    \author Hui Xue
*/

#ifndef hoImageRegHomogenousTransformation_H_
#define hoImageRegHomogenousTransformation_H_

#pragma once

#include "hoImageRegParametricTransformation.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"

namespace Gadgetron {

    /// Homogenous transformation
    template<typename ValueType, unsigned int D> 
    class hoImageRegHomogenousTransformation : public hoImageRegParametricTransformation<ValueType, D, D>
    {
    public:

        typedef hoImageRegParametricTransformation<ValueType, D, D> BaseClass;
        typedef hoImageRegHomogenousTransformation<ValueType, D> Self;

        typedef ValueType T;

        typedef typename BaseClass::input_point_type input_point_type;
        typedef typename BaseClass::output_point_type output_point_type;

        typedef typename BaseClass::jacobian_parameter_type jacobian_parameter_type;
        typedef typename BaseClass::jacobian_position_type jacobian_position_type;

        typedef typename BaseClass::ParaStatus ParaStatus;
        typedef typename BaseClass::ParaStatusType ParaStatusType;

        hoImageRegHomogenousTransformation();
        virtual ~hoImageRegHomogenousTransformation();

        // get/set the ith parameter
        virtual ValueType get_parameter(size_t i) const;
        virtual void set_parameter(size_t i, ValueType v);

        virtual bool invertTransformation();

        virtual bool setIdentity();

        virtual bool transform(const T* pt_in, T* pt_out) const;

        virtual bool transform(const T& xi, const T& yi, T& xo, T& yo) const;

        virtual bool transform(const T& xi, const T& yi, const T& zi, T& xo, T& yo, T& zo) const;

        virtual bool transform(const size_t* pt_in, T* pt_out) const;
        virtual bool transform(const size_t* pt_in, size_t N, T* pt_out) const;
        virtual bool transform(const size_t& xi, const size_t& yi, T& xo, T& yo) const;
        virtual bool transform(const size_t* xi, const size_t* yi, size_t N, T* xo, T* yo) const;
        virtual bool transform(const size_t& xi, const size_t& yi, const size_t& zi, T& xo, T& yo, T& zo) const;
        virtual bool transform(const size_t* xi, const size_t* yi, const size_t* zi, size_t N, T* xo, T* yo, T* zo) const;

        /// compute jacobian matrix to parameters
        /// D*num_parameters_ matrix
        virtual bool jacobianParameter(const input_point_type& pos, jacobian_parameter_type& jac);

        /// compute jacobian matrix to spatial position
        /// D*D matrix
        virtual bool jacobianPosition(const input_point_type& pos, jacobian_position_type& jac);

        virtual void print(std::ostream& os) const;
        virtual void printTransform(std::ostream& os) const;

        virtual std::string transformationName() const
        {
            return std::string("hoImageRegHomogenousTransformation"); 
        }

        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

    protected:

        using BaseClass::num_parameters_;
        using BaseClass::para_status_;

        /// transformation matrix
        hoMatrix<ValueType> matrix_;
    };

    template <typename ValueType, unsigned int D> 
    hoImageRegHomogenousTransformation<ValueType, D>::hoImageRegHomogenousTransformation() : BaseClass()
    {
        num_parameters_ = D*(D+1);
        para_status_.resize(num_parameters_, BaseClass::Active);

        GADGET_CHECK_THROW(matrix_.createMatrix(D+1, D+1));
        GADGET_CHECK_THROW(matrix_.setIdentity());
    }

    template <typename ValueType, unsigned int D> 
    hoImageRegHomogenousTransformation<ValueType, D>::~hoImageRegHomogenousTransformation()
    {
    }

    template <typename ValueType, unsigned int D> 
    inline ValueType hoImageRegHomogenousTransformation<ValueType, D>::get_parameter(size_t i) const
    {
        GADGET_DEBUG_CHECK_THROW(i<num_parameters_);
        return matrix_( i/(D+1), i%(D+1) );
    }

    template <typename ValueType, unsigned int D> 
    inline void hoImageRegHomogenousTransformation<ValueType, D>::set_parameter(size_t i, ValueType v)
    {
        GADGET_DEBUG_CHECK_THROW(i<num_parameters_);
        matrix_( i/(D+1), i%(D+1) ) = v;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::invertTransformation()
    {
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::invert(matrix_) );
        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::setIdentity()
    {
        GADGET_CHECK_RETURN_FALSE( matrix_.setIdentity() );
        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::transform(const T* pt_in, T* pt_out) const
    {
        try
        {
            unsigned int ii, jj;
            for ( ii=0; ii<D; ii++ )
            {
                pt_out[ii] = 0;
                for ( jj=0; jj<D; jj++ )
                {
                    pt_out[ii] += matrix_(ii, jj) * pt_in[jj];
                }
                pt_out[ii] += matrix_(ii, D);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegHomogenousTransformation<ValueType, D>::transform(const T* pt_in, T* pt_out) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::transform(const T& xi, const T& yi, T& xo, T& yo) const
    {
        try
        {
            xo = matrix_(0, 0)*xi + matrix_(0, 1)*yi + matrix_(0, 2);
            yo = matrix_(1, 0)*xi + matrix_(1, 1)*yi + matrix_(1, 2);
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegHomogenousTransformation<ValueType, D>::transform(const T& xi, const T& yi, T& xo, T& yo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::transform(const T& xi, const T& yi, const T& zi, T& xo, T& yo, T& zo) const
    {
        try
        {
            xo = matrix_(0, 0)*xi + matrix_(0, 1)*yi + matrix_(0, 2)*zi + matrix_(0, 3);
            yo = matrix_(1, 0)*xi + matrix_(1, 1)*yi + matrix_(1, 2)*zi + matrix_(1, 3);
            zo = matrix_(2, 0)*xi + matrix_(2, 1)*yi + matrix_(2, 2)*zi + matrix_(2, 3);
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegHomogenousTransformation<ValueType, D>::transform(const T& xi, const T& yi, const T& zi, T& xo, T& yo, T& zo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t* pt_in, T* pt_out) const
    {
        try
        {
            unsigned int ii, jj;
            for ( ii=0; ii<D; ii++ )
            {
                pt_out[ii] = 0;
                for ( jj=0; jj<D; jj++ )
                {
                    pt_out[ii] += matrix_(ii, jj) * pt_in[jj];
                }
                pt_out[ii] += matrix_(ii, D);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t* pt_in, T* pt_out) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t* pt_in, size_t N, T* pt_out) const
    {
        try
        {
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(pt_in, pt_out, N)
            for ( ii=0; ii<(long long)N; ii++ )
            {
                this->transform(pt_in+ii*D, pt_out+ii*D);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t* pt_in, size_t N, T* pt_out) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t& xi, const size_t& yi, T& xo, T& yo) const
    {
        try
        {
            xo = matrix_(0, 0)*xi + matrix_(0, 1)*yi + matrix_(0, 2);
            yo = matrix_(1, 0)*xi + matrix_(1, 1)*yi + matrix_(1, 2);
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t& xi, const size_t& yi, T& xo, T& yo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t* xi, const size_t* yi, size_t N, T* xo, T* yo) const
    {
        try
        {
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(xi, yi, xo, yo, N)
            for ( ii=0; ii<(long long)N; ii++ )
            {
                this->transform(xi[ii], yi[ii], xo[ii], yo[ii]);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t* xi, const size_t* yi, size_t N, T* xo, T* yo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t& xi, const size_t& yi, const size_t& zi, T& xo, T& yo, T& zo) const
    {
        try
        {
            xo = matrix_(0, 0)*xi + matrix_(0, 1)*yi + matrix_(0, 2)*zi + matrix_(0, 3);
            yo = matrix_(1, 0)*xi + matrix_(1, 1)*yi + matrix_(1, 2)*zi + matrix_(1, 3);
            zo = matrix_(2, 0)*xi + matrix_(2, 1)*yi + matrix_(2, 2)*zi + matrix_(2, 3);
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t& xi, const size_t& yi, const size_t& zi, T& xo, T& yo, T& zo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t* xi, const size_t* yi, const size_t* zi, size_t N, T* xo, T* yo, T* zo) const
    {
        try
        {
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(xi, yi, zi, xo, yo, zo, N)
            for ( ii=0; ii<(long long)N; ii++ )
            {
                this->transform(xi[ii], yi[ii], zi[ii], xo[ii], yo[ii], zo[ii]);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegHomogenousTransformation<ValueType, D>::transform(const size_t* xi, const size_t* yi, const size_t* zi, size_t N, T* xo, T* yo, T* zo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::jacobianParameter(const input_point_type& pos, jacobian_parameter_type& jac)
    {
        try
        {
            jac.createMatrix(D, num_parameters_);
            Gadgetron::clear(jac);

            if ( D == 2 )
            {
                jac(0, 0) = pos(0);
                jac(0, 1) = pos(1);
                jac(0, 2) = 1;

                jac(1, 3) = pos(0);
                jac(1, 4) = pos(1);
                jac(1, 5) = 1;
            }
            else if ( D == 3 )
            {
                jac(0, 0) = pos(0);
                jac(0, 1) = pos(1);
                jac(0, 2) = pos(2);
                jac(0, 3) = 1;

                jac(1, 4) = pos(0);
                jac(1, 5) = pos(1);
                jac(1, 6) = pos(2);
                jac(1, 7) = 1;

                jac(2, 8)  = pos(0);
                jac(2, 9)  = pos(1);
                jac(2, 10) = pos(2);
                jac(2, 11) = 1;
            }
            else
            {
                unsigned int ii, jj;
                for ( ii=0; ii<D; ii++ )
                {
                    for ( jj=0; jj<D; jj++ )
                    {
                        jac(ii, ii*(D+1)+jj) = pos(jj);
                    }

                    jac(ii, ii*(D+1)+D) = 1;
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegHomogenousTransformation<ValueType, D>::jacobianParameter(const input_point_type& pos, jacobian_parameter_type& jac) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    bool hoImageRegHomogenousTransformation<ValueType, D>::jacobianPosition(const input_point_type& pos, jacobian_position_type& jac)
    {
        try
        {
            jac.createMatrix(D, D);
            Gadgetron::clear(jac);

            if ( D == 2 )
            {
                jac(0, 0) = matrix_(0, 0);
                jac(0, 1) = matrix_(0, 1);
                jac(1, 0) = matrix_(1, 0);
                jac(1, 1) = matrix_(1, 1);
            }
            else if ( D == 3 )
            {
                jac(0, 0) = matrix_(0, 0);
                jac(0, 1) = matrix_(0, 1);
                jac(0, 2) = matrix_(0, 2);

                jac(1, 0) = matrix_(1, 0);
                jac(1, 1) = matrix_(1, 1);
                jac(1, 2) = matrix_(1, 2);

                jac(2, 0) = matrix_(2, 0);
                jac(2, 1) = matrix_(2, 1);
                jac(2, 2) = matrix_(2, 2);
            }
            else
            {
                unsigned int ii, jj;
                for ( ii=0; ii<D; ii++ )
                {
                    for ( jj=0; jj<D; jj++ )
                    {
                        jac(ii, jj) = matrix_(ii, jj);
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegHomogenousTransformation<ValueType, D>::jacobianPosition(const input_point_type& pos, jacobian_position_type& jac) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D> 
    void hoImageRegHomogenousTransformation<ValueType, D>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron homogenous transformation -------------" << endl;
        os << "Input dimension is : " << D << endl;
        os << "Output dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(T).name());
        os << "Transformation data type is : " << elemTypeName << std::endl;
        os << "Number of parameters is : " << num_parameters_ << endl;

        size_t i;
        os << "Status of parameters: " << endl;
        for ( i=0; i<this->num_parameters_; i++ )
        {
            os << "Para " << i << " : \t";
            if ( para_status_[i] == BaseClass::Active )
            {
                os << "Active";
            }
            else if ( para_status_[i] == BaseClass::Inactive )
            {
                os << "Inactive";
            }
            else
            {
                os << "Unknown";
            }
            os << endl;
        }

        os << "Transformation: " << endl;
        this->printTransform(os);
    }

    template <typename ValueType, unsigned int D> 
    void hoImageRegHomogenousTransformation<ValueType, D>::printTransform(std::ostream& os) const
    {
        using namespace std;

        size_t i;
        os << "[ ";
        for ( i=0; i<this->num_parameters_; i++ )
        {
            os << this->get_parameter(i) << " \t";
        }
        os << " ]" << endl;
    }
}
#endif // hoImageRegHomogenousTransformation_H_
