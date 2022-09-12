/** \file   hoImageRegRigid3DTransformation.h
    \brief  Define the class for the rigid 2D transformation in gadgetron registration
            Three parameters are translation along x and y and roation along z (tx, ty, rz)
    \author Hui Xue
*/

#ifndef hoImageRegRigid3DTransformation_H_
#define hoImageRegRigid3DTransformation_H_

#pragma once

#include "hoImageRegHomogenousTransformation.h"
#include <cmath>

namespace Gadgetron {

    /// Homogenous transformation
    template<typename ValueType> 
    class hoImageRegRigid3DTransformation : public hoImageRegHomogenousTransformation<ValueType, 3>
    {
    public:

        typedef hoImageRegParametricTransformation<ValueType, 3, 3> ParaTransformBaseClass;
        typedef hoImageRegHomogenousTransformation<ValueType, 3> BaseClass;
        typedef hoImageRegRigid3DTransformation<ValueType> Self;

        typedef ValueType T;

        typedef typename BaseClass::input_point_type input_point_type;
        typedef typename BaseClass::output_point_type output_point_type;

        typedef typename BaseClass::jacobian_parameter_type jacobian_parameter_type;
        typedef typename BaseClass::jacobian_position_type jacobian_position_type;

        typedef typename BaseClass::ParaStatus ParaStatus;
        typedef typename BaseClass::ParaStatusType ParaStatusType;

        hoImageRegRigid3DTransformation();
        virtual ~hoImageRegRigid3DTransformation();

        // get/set the ith parameter
        virtual ValueType get_parameter(size_t i) const;
        virtual void set_parameter(size_t i, ValueType v);

        virtual bool invertTransformation();

        virtual bool setIdentity();

        // get/set the translation and rotation
        ValueType get_tx() const;
        ValueType get_ty() const;
        ValueType get_tz() const;
        ValueType get_rx() const;
        ValueType get_ry() const;
        ValueType get_rz() const;

        void set_tx(ValueType tx);
        void set_ty(ValueType ty);
        void set_tz(ValueType tz);
        void set_rx(ValueType rx);
        void set_ry(ValueType ry);
        void set_rz(ValueType rz);

        void set_tx_ty_tz(ValueType tx, ValueType ty, ValueType tz);
        void set_rx_ry_rz(ValueType rx, ValueType ry, ValueType rz);

        /// compute the transformation matrix
        bool updateTransformationMatrix(ValueType tx, ValueType ty, ValueType tz, ValueType rx, ValueType ry, ValueType rz, hoMatrix<T>& matrix);
        bool extractParametersFromTransformationMatrix(const hoMatrix<T>& matrix, ValueType& tx, ValueType& ty, ValueType& tz, ValueType& rx, ValueType& ry, ValueType& rz);

        virtual bool adjustForResolutionPyramid(const hoMatrix<ValueType>& sourceI2W, const hoMatrix<ValueType>& targetI2W);

        /// compute jacobian matrix to parameters
        /// D*num_parameters_ matrix
        virtual bool jacobianParameter(const input_point_type& pos, jacobian_parameter_type& jac);

        virtual void print(std::ostream& os) const;
        virtual void printTransform(std::ostream& os) const;

        virtual std::string transformationName() const
        {
            return std::string("hoImageRegRigid3DTransformation"); 
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

        /// translation along x, y and z
        ValueType tx_;
        ValueType ty_;
        ValueType tz_;
        /// rotation along x, y, and z, in degree
        ValueType rx_;
        ValueType ry_;
        ValueType rz_;
    };

    template <typename ValueType> 
    hoImageRegRigid3DTransformation<ValueType>::hoImageRegRigid3DTransformation() : BaseClass()
    {
        num_parameters_ = 6;
        para_status_.resize(num_parameters_, BaseClass::Active);

        GADGET_CHECK_THROW(matrix_.createMatrix(4, 4));
        GADGET_CHECK_THROW(matrix_.setIdentity());

        tx_ = 0;
        ty_ = 0;
        tz_ = 0;
        rx_ = 0;
        ry_ = 0;
        rz_ = 0;
    }

    template <typename ValueType> 
    hoImageRegRigid3DTransformation<ValueType>::~hoImageRegRigid3DTransformation()
    {
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid3DTransformation<ValueType>::get_parameter(size_t i) const
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
        else if ( i == 2 )
        {
            return tz_;
        }
        else if ( i == 3 )
        {
            return rx_;
        }
        else if ( i == 4 )
        {
            return ry_;
        }
        else if ( i == 5 )
        {
            return rz_;
        }

        return 0;
    }

    template <typename ValueType> 
    inline void hoImageRegRigid3DTransformation<ValueType>::set_parameter(size_t i, ValueType v)
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
        else if ( i == 2 )
        {
            tz_ = v;
        }
        else if ( i == 3 )
        {
            rx_ = v;
        }
        else if ( i == 4 )
        {
            ry_ = v;
        }
        else if ( i == 5 )
        {
            rz_ = v;
        }

        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline bool hoImageRegRigid3DTransformation<ValueType>::invertTransformation()
    {
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::invert(matrix_) );
        GADGET_CHECK_RETURN_FALSE( this->extractParametersFromTransformationMatrix(matrix_, tx_, ty_, tz_, rx_, ry_, rz_) );
        return true;
    }

    template <typename ValueType> 
    inline bool hoImageRegRigid3DTransformation<ValueType>::setIdentity()
    {
        GADGET_CHECK_RETURN_FALSE( matrix_.setIdentity() );
        tx_ = 0;
        ty_ = 0;
        tz_ = 0;
        rx_ = 0;
        ry_ = 0;
        rz_ = 0;
        return true;
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid3DTransformation<ValueType>::get_tx() const
    {
        return tx_;
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid3DTransformation<ValueType>::get_ty() const
    {
        return ty_;
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid3DTransformation<ValueType>::get_tz() const
    {
        return tz_;
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid3DTransformation<ValueType>::get_rx() const
    {
        return rx_;
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid3DTransformation<ValueType>::get_ry() const
    {
        return ry_;
    }

    template <typename ValueType> 
    inline ValueType hoImageRegRigid3DTransformation<ValueType>::get_rz() const
    {
        return rz_;
    }

    template <typename ValueType> 
    inline void hoImageRegRigid3DTransformation<ValueType>::set_tx(ValueType tx)
    {
        tx_ = tx;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid3DTransformation<ValueType>::set_ty(ValueType ty)
    {
        ty_ = ty;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid3DTransformation<ValueType>::set_tz(ValueType tz)
    {
        tz_ = tz;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid3DTransformation<ValueType>::set_rx(ValueType rx)
    {
        rx_ = rx;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid3DTransformation<ValueType>::set_ry(ValueType ry)
    {
        ry_ = ry;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid3DTransformation<ValueType>::set_rz(ValueType rz)
    {
        rz_ = rz;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid3DTransformation<ValueType>::set_tx_ty_tz(ValueType tx, ValueType ty, ValueType tz)
    {
        tx_ = tx;
        ty_ = ty;
        tz_ = tz;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
    }

    template <typename ValueType> 
    inline void hoImageRegRigid3DTransformation<ValueType>::set_rx_ry_rz(ValueType rx, ValueType ry, ValueType rz)
    {
        rx_ = rx;
        ry_ = ry;
        rz_ = rz;
        GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
    }

    template <typename ValueType> 
    bool hoImageRegRigid3DTransformation<ValueType>::updateTransformationMatrix(ValueType tx, ValueType ty, ValueType tz, ValueType rx, ValueType ry, ValueType rz, hoMatrix<T>& matrix)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE( matrix.createMatrix(4, 4) );

            double cosrx = std::cos(rx*M_PI/180.0);
            double sinrx = std::sin(rx*M_PI/180.0);

            double cosry = std::cos(ry*M_PI/180.0);
            double sinry = std::sin(ry*M_PI/180.0);

            double cosrz = std::cos(rz*M_PI/180.0);
            double sinrz = std::sin(rz*M_PI/180.0);

            matrix(0, 0) = cosry*cosrz;                         matrix(0, 1) = cosry*sinrz;                             matrix(0, 2) = -sinry;           matrix(0, 3) = tx;
            matrix(1, 0) = sinrx*sinry*cosrz-cosrx*sinrz;       matrix(1, 1) = sinrx*sinry*sinrz+cosrx*cosrz;           matrix(1, 2) = sinrx*cosry;      matrix(1, 3) = ty;
            matrix(2, 0) = cosrx*sinry*cosrz+sinrx*sinrz;       matrix(2, 1) = cosrx*sinry*sinrz-sinrx*cosrz;           matrix(2, 2) = cosrx*cosry;      matrix(2, 3) = tz;
            matrix(3, 0) = 0;                                   matrix(3, 1) = 0;                                       matrix(3, 2) = 0;                matrix(3, 3) = 1;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegRigid3DTransformation<ValueType>::updateTransformationMatrix(ValueType tx, ValueType ty, ValueType rz, hoMatrix<T>& matrix) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType> 
    bool hoImageRegRigid3DTransformation<ValueType>::extractParametersFromTransformationMatrix(const hoMatrix<T>& matrix, ValueType& tx, ValueType& ty, ValueType& tz, ValueType& rx, ValueType& ry, ValueType& rz)
    {
        try
        {
            ry_ = asin(-1 * matrix_(0, 2));

            if ( std::abs( std::cos(ry_) ) > 1e-6 )
            {
                rx_ = atan2(matrix_(1, 2), matrix_(2, 2));
                rz_ = atan2(matrix_(0, 1), matrix_(0, 0));
            } 
            else 
            { 
                rx_ = atan2(-1.0*matrix_(0, 2)*matrix_(1, 0), -1.0*matrix_(0, 2)*matrix_(2, 0)); 
                rz_ = 0;
            }

            tx_ = matrix_(0, 3);
            ty_ = matrix_(1, 3);
            tz_ = matrix_(2, 3);
            rx_ *= 180.0/M_PI;
            ry_ *= 180.0/M_PI;
            rz_ *= 180.0/M_PI;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegRigid3DTransformation<ValueType>::extractParametersFromTransformationMatrix(const hoMatrix<T>& matrix, ValueType& tx, ValueType& ty, ValueType& tz, ValueType& rx, ValueType& ry, ValueType& rz) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType> 
    bool hoImageRegRigid3DTransformation<ValueType>::jacobianParameter(const input_point_type& pos, jacobian_parameter_type& jac)
    {
        try
        {
            jac.createMatrix(3, num_parameters_);
            Gadgetron::clear(jac);

            double cosrx = std::cos(rx_*M_PI/180.0);
            double sinrx = std::sin(rx_*M_PI/180.0);

            double cosry = std::cos(ry_*M_PI/180.0);
            double sinry = std::sin(ry_*M_PI/180.0);

            double cosrz = std::cos(rz_*M_PI/180.0);
            double sinrz = std::sin(rz_*M_PI/180.0);

            jac(0, 0) = 1;
            jac(0, 1) = 0;
            jac(0, 2) = 0;
            jac(0, 3) = 0;
            jac(0, 4) = -sinry*cosrz*pos(0)-sinry*sinrz*pos(1)-cosry*pos(2);
            jac(0, 5) = -cosry*sinrz*pos(0)+cosry*cosrz*pos(1);

            jac(1, 0) = 0;
            jac(1, 1) = 1;
            jac(1, 2) = 0;
            jac(1, 3) = (cosrx*sinry*cosrz+sinrx*sinrz) *pos(0)             + (cosrx*sinry*sinrz-sinrx*cosrz)   *pos(1)        + cosrx*cosry*pos(2);
            jac(1, 4) = (sinrx*cosry*cosrz)             *pos(0)             + (sinrx*cosry*sinrz)               *pos(1)        - sinrx*sinry*pos(2);
            jac(1, 5) = (-sinrx*sinry*sinrz-cosrx*cosrz)*pos(0)             + (sinrx*sinry*cosrz-cosrx*sinrz)   *pos(1);

            jac(2, 0) = 0;
            jac(2, 1) = 0;
            jac(2, 2) = 1;
            jac(2, 3) = (-sinrx*sinry*cosrz+cosrx*sinrz)*pos(0)             + (-sinrx*sinry*sinrz-cosrx*cosrz)  *pos(1)         - sinrx*cosry*pos(2);
            jac(2, 4) = (cosrx*cosry*cosrz)             *pos(0)             + (cosrx*cosry*sinrz)               *pos(1)         - cosrx*sinry*pos(2);
            jac(2, 5) = (cosrx*sinry*-sinrz+sinrx*cosrz)*pos(0)             + (cosrx*sinry*cosrz+sinrx*sinrz)   *pos(1);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happen in hoImageRegRigid3DTransformation<ValueType>::jacobianParameter(const input_point_type& pos, jacobian_parameter_type& jac) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType> 
    bool hoImageRegRigid3DTransformation<ValueType>::adjustForResolutionPyramid(const hoMatrix<ValueType>& sourceI2W, const hoMatrix<ValueType>& targetI2W)
    {
        try
        {
            hoNDImage<ValueType, 3> source;
            source.set_image_to_world_matrix(sourceI2W);

            hoNDImage<ValueType, 3> target;
            target.set_image_to_world_matrix(targetI2W);

            tx_ *= source.get_pixel_size(0)/target.get_pixel_size(0);
            ty_ *= source.get_pixel_size(1)/target.get_pixel_size(1);
            tz_ *= source.get_pixel_size(2)/target.get_pixel_size(2);

            GADGET_CHECK_THROW(this->updateTransformationMatrix(tx_, ty_, tz_, rx_, ry_, rz_, matrix_));
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegRigid3DTransformation<ValueType>::adjustForResolutionPyramid(const hoMatrix<ValueType>& sourceI2W, const hoMatrix<ValueType>& targetI2W) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType> 
    void hoImageRegRigid3DTransformation<ValueType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron rigid 3D transformation -------------" << endl;
        std::string elemTypeName = std::string(typeid(T).name());
        os << "Transformation data type is : " << elemTypeName << std::endl;
        os << "Number of parameters is : " << num_parameters_ << endl;

        size_t i;
        os << "Status of parameters [tx ty tz rx ry rz] : " << endl;
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
    void hoImageRegRigid3DTransformation<ValueType>::printTransform(std::ostream& os) const
    {
        using namespace std;

        size_t i;
        os << "[tx ty tz rx ry rz] = [ ";
        for ( i=0; i<this->num_parameters_; i++ )
        {
            os << this->get_parameter(i) << " \t";
        }
        os << " ]" << endl;
    }
}
#endif // hoImageRegRigid2DTransformation_H_
