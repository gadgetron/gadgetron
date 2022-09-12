/** \file   hoImageRegRigid2DTransformation.h
    \brief  Define the class for the rigid 2D transformation in gadgetron registration
            Three parameters are translation along x and y and roation along z (tx, ty, rz)
    \author Hui Xue
*/

#ifndef hoImageRegRigid2DTransformation_H_
#define hoImageRegRigid2DTransformation_H_

#pragma once

#include "hoImageRegHomogenousTransformation.h"
#include <cmath>

namespace Gadgetron {

    /// Homogenous transformation
    template<typename ValueType> 
    class hoImageRegRigid2DTransformation : public hoImageRegHomogenousTransformation<ValueType, 2>
    {
    public:

        typedef hoImageRegParametricTransformation<ValueType, 2, 2> ParaTransformBaseClass;
        typedef hoImageRegHomogenousTransformation<ValueType, 2> BaseClass;
        typedef hoImageRegRigid2DTransformation<ValueType> Self;

        typedef ValueType T;

        typedef typename BaseClass::input_point_type input_point_type;
        typedef typename BaseClass::output_point_type output_point_type;

        typedef typename BaseClass::jacobian_parameter_type jacobian_parameter_type;
        typedef typename BaseClass::jacobian_position_type jacobian_position_type;

        typedef typename BaseClass::ParaStatus ParaStatus;
        typedef typename BaseClass::ParaStatusType ParaStatusType;

        hoImageRegRigid2DTransformation();
        virtual ~hoImageRegRigid2DTransformation();

        // get/set the ith parameter
        virtual ValueType get_parameter(size_t i) const;
        virtual void set_parameter(size_t i, ValueType v);

        virtual bool invertTransformation();

        virtual bool setIdentity();

        // get/set the translation and rotation
        ValueType get_tx() const;
        ValueType get_ty() const;
        ValueType get_rz() const;

        void set_tx(ValueType tx);
        void set_ty(ValueType ty);
        void set_rz(ValueType rz);

        void set_tx_ty(ValueType tx, ValueType ty);
        void set_tx_ty_rz(ValueType tx, ValueType ty, ValueType rz);

        /// compute the transformation matrix
        bool updateTransformationMatrix(ValueType tx, ValueType ty, ValueType rz, hoMatrix<T>& matrix);
        bool extractParametersFromTransformationMatrix(const hoMatrix<T>& matrix, ValueType& tx, ValueType& ty, ValueType& rz);

        virtual bool adjustForResolutionPyramid(const hoMatrix<ValueType>& sourceI2W, const hoMatrix<ValueType>& targetI2W);

        /// compute jacobian matrix to parameters
        /// D*num_parameters_ matrix
        virtual bool jacobianParameter(const input_point_type& pos, jacobian_parameter_type& jac);

        virtual void print(std::ostream& os) const;
        virtual void printTransform(std::ostream& os) const;

        virtual std::string transformationName() const
        {
            return std::string("hoImageRegRigid2DTransformation"); 
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
        using BaseClass::matrix_;

        /// translation along x and y
        ValueType tx_;
        ValueType ty_;
        /// rotation along z, in degree
        ValueType rz_;
    };

    template <typename ValueType> 
    hoImageRegRigid2DTransformation<ValueType>::hoImageRegRigid2DTransformation() : BaseClass()
    {
        num_parameters_ = 3;
        para_status_.resize(num_parameters_, ParaTransformBaseClass::Active);

        GADGET_CHECK_THROW(matrix_.createMatrix(3, 3));
        GADGET_CHECK_THROW(matrix_.setIdentity());

        tx_ = 0;
        ty_ = 0;
        rz_ = 0;
    }

    template <typename ValueType> 
    hoImageRegRigid2DTransformation<ValueType>::~hoImageRegRigid2DTransformation()
    {
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid2DTransformation<ValueType>::get_parameter(size_t i) const
    {
        GADGET_DEBUG_CHECK_THROW(i<num_parameters_);
        if ( i == 0 )
        {
            return tx_;
        }
        else if ( i == 1 )
        {
            return ty_;
        }
        else
        {
            return rz_;
        }
    }

    template <typename ValueType> 
    inline void hoImageRegRigid2DTransformation<ValueType>::set_parameter(size_t i, ValueType v)
    {
        GADGET_DEBUG_CHECK_THROW(i<num_parameters_);
        if ( i == 0 )
        {
            tx_ = v;
        }
        else if ( i == 1 )
        {
            ty_ = v;
        }
        else
        {
            rz_ = v;
        }

        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline bool hoImageRegRigid2DTransformation<ValueType>::invertTransformation()
    {
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::invert(matrix_) );
        GADGET_CHECK_RETURN_FALSE( this->extractParametersFromTransformationMatrix(matrix_, tx_, ty_, rz_) );
        return true;
    }

    template <typename ValueType> 
    inline bool hoImageRegRigid2DTransformation<ValueType>::setIdentity()
    {
        GADGET_CHECK_RETURN_FALSE( matrix_.setIdentity() );
        tx_ = 0;
        ty_ = 0;
        rz_ = 0;
        return true;
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid2DTransformation<ValueType>::get_tx() const
    {
        return tx_;
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid2DTransformation<ValueType>::get_ty() const
    {
        return ty_;
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid2DTransformation<ValueType>::get_rz() const
    {
        return rz_;
    }

    template <typename ValueType> 
    inline void hoImageRegRigid2DTransformation<ValueType>::set_tx(ValueType tx)
    {
        tx_ = tx;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid2DTransformation<ValueType>::set_ty(ValueType ty)
    {
        ty_ = ty;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid2DTransformation<ValueType>::set_rz(ValueType rz)
    {
        rz_ = rz;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid2DTransformation<ValueType>::set_tx_ty(ValueType tx, ValueType ty)
    {
        tx_ = tx;
        ty_ = ty;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid2DTransformation<ValueType>::set_tx_ty_rz(ValueType tx, ValueType ty, ValueType rz)
    {
        tx_ = tx;
        ty_ = ty;
        rz_ = rz;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, rz_, matrix_));
    }

    template <typename ValueType> 
    bool hoImageRegRigid2DTransformation<ValueType>::updateTransformationMatrix(ValueType tx, ValueType ty, ValueType rz, hoMatrix<T>& matrix)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE( matrix.createMatrix(3, 3) );

            ValueType cosrz = std::cos(rz*M_PI/180.0);
            ValueType sinrz = std::sin(rz*M_PI/180.0);

            matrix(0, 0) = cosrz;  matrix(0, 1) = sinrz; matrix(0, 2) = tx;
            matrix(1, 0) = -sinrz; matrix(1, 1) = cosrz; matrix(1, 2) = ty;
            matrix(2, 0) = 0;      matrix(2, 1) = 0;     matrix(2, 2) = 1;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegRigid2DTransformation<ValueType>::updateTransformationMatrix(ValueType tx, ValueType ty, ValueType rz, hoMatrix<T>& matrix) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType> 
    bool hoImageRegRigid2DTransformation<ValueType>::extractParametersFromTransformationMatrix(const hoMatrix<T>& matrix, ValueType& tx, ValueType& ty, ValueType& rz)
    {
        try
        {
            double cosrz = matrix(0, 0);
            double sinrz = matrix(0, 1);

            if ( cosrz >= 0 ) // rz is [-PI/2 PI/2]
            {
                rz = std::asin(sinrz);
            }
            else
            {
                rz = std::acos(cosrz);
                if ( sinrz < 0) rz *= -1; // [-PI -PI/2]
            }

            tx = matrix(0, 2);
            ty = matrix(1, 2);
            rz *= 180.0/M_PI;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegRigid2DTransformation<ValueType>::extractParametersFromTransformationMatrix(const hoMatrix<T>& matrix, ValueType& tx, ValueType& ty, ValueType& rz) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType> 
    bool hoImageRegRigid2DTransformation<ValueType>::jacobianParameter(const input_point_type& pos, jacobian_parameter_type& jac)
    {
        try
        {
            jac.createMatrix(2, num_parameters_);
            Gadgetron::clear(jac);

            double cosrz = matrix_(0, 0);
            double sinrz = matrix_(0, 1);

            jac(0, 0) = 1;
            jac(0, 1) = 0;
            jac(0, 2) = -sinrz*pos(0) + cosrz*pos(1);

            jac(1, 0) = 0;
            jac(1, 1) = 1;
            jac(1, 2) = -cosrz*pos(0) - sinrz*pos(1);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegRigid2DTransformation<ValueType>::jacobianParameter(const input_point_type& pos, jacobian_parameter_type& jac) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType> 
    bool hoImageRegRigid2DTransformation<ValueType>::adjustForResolutionPyramid(const hoMatrix<ValueType>& sourceI2W, const hoMatrix<ValueType>& targetI2W)
    {
        try
        {
            hoNDImage<ValueType, 2> source;
            source.set_image_to_world_matrix(sourceI2W);

            hoNDImage<ValueType, 2> target;
            target.set_image_to_world_matrix(targetI2W);

            tx_ *= source.get_pixel_size(0)/target.get_pixel_size(0);
            ty_ *= source.get_pixel_size(1)/target.get_pixel_size(1);

            GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, rz_, matrix_));
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegRigid2DTransformation<ValueType>::adjustForResolutionPyramid(const hoMatrix<ValueType>& sourceI2W, const hoMatrix<ValueType>& targetI2W) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType> 
    void hoImageRegRigid2DTransformation<ValueType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron rigid 2D transformation -------------" << endl;
        std::string elemTypeName = std::string(typeid(T).name());
        os << "Transformation data type is : " << elemTypeName << std::endl;
        os << "Number of parameters is : " << num_parameters_ << endl;

        size_t i;
        os << "Status of parameters [tx ty rz] : " << endl;
        for ( i=0; i<this->num_parameters_; i++ )
        {
            os << "Para " << i << " : \t";
            if ( para_status_[i] == ParaTransformBaseClass::Active )
            {
                os << "Active";
            }
            else if ( para_status_[i] == ParaTransformBaseClass::Inactive )
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

    template <typename ValueType> 
    void hoImageRegRigid2DTransformation<ValueType>::printTransform(std::ostream& os) const
    {
        using namespace std;

        size_t i;
        os << "[tx ty rz] = [ ";
        for ( i=0; i<this->num_parameters_; i++ )
        {
            os << this->get_parameter(i) << " \t";
        }
        os << " ]" << endl;
    }
}
#endif // hoImageRegRigid2DTransformation_H_
