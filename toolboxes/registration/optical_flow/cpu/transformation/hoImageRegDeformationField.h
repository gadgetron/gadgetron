/** \file   hoImageRegDeformationField.h
    \brief  Define the geometry transformation using deformation filed
    \author Hui Xue
*/

#ifndef hoImageRegDeformationField_H_
#define hoImageRegDeformationField_H_

#pragma once

#include "hoImageRegNonParametricTransformation.h"

namespace Gadgetron {

    /// deformation field is defined as hoNDImage
    /// the deformation field can be accessed on image pixels
    /// if the non-integer image pixels are used to access deformaiton field, an image interpolator is used
    /// linear interpolator is used for deformation field
    /// the unit of stored deformation field is in pixel, not in world coordinates
    template<typename ValueType, unsigned int D>
    class  hoImageRegDeformationField: public hoImageRegNonParametricTransformation<ValueType, D, D>
    {
    public:

        typedef hoImageRegTransformation<ValueType, D, D> Self;
        typedef hoImageRegNonParametricTransformation<ValueType, D, D> BaseClass;

        typedef ValueType T;

        typedef typename BaseClass::input_point_type input_point_type;
        typedef typename BaseClass::output_point_type output_point_type;

        typedef typename BaseClass::jacobian_position_type jacobian_position_type;

        typedef hoNDImage<T, D> DeformationFieldType;

        typedef typename DeformationFieldType::coord_type coord_type;

        typedef typename DeformationFieldType::axis_type axis_type;

        typedef hoNDBoundaryHandler<DeformationFieldType> BoundHanlderType;
        typedef hoNDInterpolator<DeformationFieldType> InterpolatorType;

        typedef hoNDInterpolatorLinear<DeformationFieldType> DefaultInterpolatorType;
        typedef hoNDBoundaryHandlerBorderValue<DeformationFieldType> DefaultBoundHanlderType;

        hoImageRegDeformationField();
        hoImageRegDeformationField(const std::vector<size_t>& dimensions);
        hoImageRegDeformationField(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis);
        hoImageRegDeformationField(const hoNDImage<ValueType, D>& im);

        virtual ~hoImageRegDeformationField();

        virtual bool invertTransformation();

        virtual bool setIdentity();

        /// update the internal status after the deformation fields are changed
        virtual bool update();

        /// transform a point
        /// the point is in the non-integer image pixel indexes
        /// image interpolator is used
        virtual bool transform(const T* pt_in, T* pt_out) const;
        virtual bool transform(const T& xi, const T& yi, T& xo, T& yo) const;
        virtual bool transform(const T& xi, const T& yi, const T& zi, T& xo, T& yo, T& zo) const;

        /// transform a point
        /// the point is in the integer image pixel indexes
        /// image interpolator is not used
        /// pt_in, pt_out stores a point as an array
        virtual bool transform(const size_t* pt_in, T* pt_out) const;
        virtual bool transform(const size_t* pt_in, size_t N, T* pt_out) const;

        /// for 2D - 2D transformation
        virtual bool transform(const size_t& xi, const size_t& yi, T& xo, T& yo) const;
        virtual bool transform(const size_t* xi, const size_t* yi, size_t N, T* xo, T* yo) const;

        /// for 3D - 3D transformation
        virtual bool transform(const size_t& xi, const size_t& yi, const size_t& zi, T& xo, T& yo, T& zo) const;
        virtual bool transform(const size_t* xi, const size_t* yi, const size_t* zi, size_t N, T* xo, T* yo, T* zo) const;

        /// compute jacobian matrix to spatial position
        /// the jacobian matrix is computed with the compensation for non-isotropic pixel sizes
        /// e.g. dxdy = ( dx(x,y+dh)*sx - dx(x, y-dh)*sx ) / (2*dh*sy); sx, sy: pixel sizes for x and y directions
        /// DOut*DIn matrix
        virtual bool jacobianPosition(const input_point_type& /*pos*/, jacobian_position_type& jac);

        /// compute jacobian matrix on the deformation grid
        /// jac is [DOut Din dimensions] array, storing the jacobian matrix for every point in the deformation field
        virtual bool jacobianPosition(hoNDArray<T>& jac, DeformationFieldType* deform_field[D], unsigned int borderWidth=1);
        virtual bool jacobianPosition(hoNDArray<T>& jac, unsigned int borderWidth=1);

        /// compute some parameters from deformation field and jacobian matrix
        /// in the world coordinate
        virtual bool analyzeJacobianAndDeformation(const hoNDArray<T>& jac, DeformationFieldType* deform_field[D], T& meanDeform, T& maxDeform, T& meanLogJac, T& maxLogJac, unsigned int borderWidth=1);
        virtual bool analyzeJacobianAndDeformation(const hoNDArray<T>& jac, T& meanDeform, T& maxDeform, T& meanLogJac, T& maxLogJac, unsigned int borderWidth=1);

        /// get/set the deformation vector on the deformation grid (image coordinate)
        /// given the index idx[DIn], output the deformation value for outDim
        T& operator()( size_t idx[D], size_t outDim );
        const T& operator()( size_t idx[D], size_t outDim ) const;

        void get(size_t idx[D], T deform[D]);
        void get(size_t x, size_t y, T& dx, T& dy);
        void get(size_t x, size_t y, size_t z, T& dx, T& dy, T& dz);

        void set(size_t idx[D], T deform[D]);
        void set(size_t x, size_t y, T dx, T dy);
        void set(size_t x, size_t y, size_t z, T dx, T dy, T dz);

        /// get/set the deformation vector on the world coordinate
        /// given the position pos[DIn], output the deformation value for outDim
        T operator()( coord_type pos[D], size_t outDim );

        void get(coord_type pos[D], T deform[D]);
        void get(coord_type px, coord_type py, T& dx, T& dy);
        void get(coord_type px, coord_type py, coord_type pz, T& dx, T& dy, T& dz);

        /// get/set interpolator
        //void getInterpolator(InterpolatorType*& interp, size_t outDim);
        //void setInterpolator(InterpolatorType* interp, size_t outDim);

        /// get/set deformation field
        void getDeformationField(DeformationFieldType*& deform, size_t outDim);
        DeformationFieldType& getDeformationField(size_t outDim) { GADGET_DEBUG_CHECK_THROW(outDim<=D); return this->deform_field_[outDim]; }

        void setDeformationField(const DeformationFieldType& deform, size_t outDim);

        /// serialize/deserialize the transformation
        virtual bool serialize(char*& buf, size_t& len) const ;
        virtual bool deserialize(char* buf, size_t& len);

        virtual void print(std::ostream& os) const;

        virtual std::string transformationName() const;

        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

    protected:

        DeformationFieldType deform_field_[D];

        //InterpolatorType* interp_[D];

        DefaultInterpolatorType* interp_default_[D];
        DefaultBoundHanlderType* bh_default_[D];
    };

    template <typename ValueType, unsigned int D>
    hoImageRegDeformationField<ValueType, D>::hoImageRegDeformationField() : BaseClass()
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            //interp_[ii] = NULL;
            bh_default_[ii] = new DefaultBoundHanlderType(deform_field_[ii]);
            interp_default_[ii] = new DefaultInterpolatorType(deform_field_[ii], *(bh_default_[ii]));
        }
    }

    template <typename ValueType, unsigned int D>
    hoImageRegDeformationField<ValueType, D>::
    hoImageRegDeformationField(const std::vector<size_t>& dimensions) : BaseClass()
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            deform_field_[ii].create(dimensions);
            memset(deform_field_[ii].get_data_ptr(), 0, deform_field_[ii].get_number_of_elements()*sizeof(T));

            //interp_[ii] = NULL;

            bh_default_[ii] = new DefaultBoundHanlderType(deform_field_[ii]);
            interp_default_[ii] = new DefaultInterpolatorType(deform_field_[ii], *(bh_default_[ii]));
        }
    }

    template <typename ValueType, unsigned int D>
    hoImageRegDeformationField<ValueType, D>::
    hoImageRegDeformationField(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis) : BaseClass()
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            deform_field_[ii].create(dimensions, pixelSize, origin, axis);
            memset(deform_field_[ii].get_data_ptr(), 0, deform_field_[ii].get_number_of_elements()*sizeof(T));

            //interp_[ii] = NULL;

            bh_default_[ii] = new DefaultBoundHanlderType(deform_field_[ii]);
            interp_default_[ii] = new DefaultInterpolatorType(deform_field_[ii], *(bh_default_[ii]));
        }
    }

    template <typename ValueType, unsigned int D>
    hoImageRegDeformationField<ValueType, D>::hoImageRegDeformationField(const hoNDImage<ValueType, D>& im) : BaseClass()
    {
        std::vector<size_t> dim;
        im.get_dimensions(dim);

        std::vector<coord_type> pixelSize;
        im.get_pixel_size(pixelSize);

        std::vector<coord_type> origin;
        im.get_origin(origin);

        axis_type axis;
        im.get_axis(axis);

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            deform_field_[ii].create(dim, pixelSize, origin, axis);
            memset(deform_field_[ii].get_data_ptr(), 0, deform_field_[ii].get_number_of_elements()*sizeof(T));

            //interp_[ii] = NULL;

            bh_default_[ii] = new DefaultBoundHanlderType(deform_field_[ii]);
            interp_default_[ii] = new DefaultInterpolatorType(deform_field_[ii], *(bh_default_[ii]));
        }
    }

    template <typename ValueType, unsigned int D>
    hoImageRegDeformationField<ValueType, D>::
    ~hoImageRegDeformationField()
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            delete bh_default_[ii];
            delete interp_default_[ii];
        }
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::invertTransformation()
    {
        /// to be implemented ...
        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::setIdentity()
    {
        try
        {
            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                memset(deform_field_[ii].get_data_ptr(), 0, deform_field_[ii].get_number_of_elements()*sizeof(T));
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::setIdentity() ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::update()
    {
        try
        {
            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                interp_default_[ii]->setArray(deform_field_[ii]);
                interp_default_[ii]->setBoundaryHandler(*bh_default_[ii]);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::update() ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::transform(const T* pt_in, T* pt_out) const
    {
        try
        {
            std::vector<coord_type> pos(D);

            int ii;
            for ( ii=0; ii<(int)D; ii++ )
            {
                pos[ii] = pt_in[ii];
            }

            #pragma omp parallel for default(none) private(ii) shared(pos, pt_out)
            for ( ii=0; ii<(int)D; ii++ )
            {
                pt_out[ii] += this->interp_default_[ii]->operator()(pos);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::transform(T* pt_in, T* pt_out) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::transform(const T& xi, const T& yi, T& xo, T& yo) const
    {
        try
        {
            xo = xi + (*interp_default_[0])(xi, yi);
            yo = yi + (*interp_default_[1])(xi, yi);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::transform(const T& xi, const T& yi, T& xo, T& yo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::transform(const T& xi, const T& yi, const T& zi, T& xo, T& yo, T& zo) const
    {
        try
        {
            xo = xi + (*interp_default_[0])(xi, yi, zi);
            yo = yi + (*interp_default_[1])(xi, yi, zi);
            zo = zi + (*interp_default_[2])(xi, yi, zi);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::transform(const T& xi, const T& yi, const T& zi, T& xo, T& yo, T& zo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::transform(const size_t* pt_in, T* pt_out) const
    {
        try
        {
            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                pt_out[ii] = pt_in[ii] + this->deform_field_[ii](pt_in);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::transform(size_t* pt_in, T* pt_out) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::transform(const size_t* pt_in, size_t N, T* pt_out) const
    {
        try
        {
            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pt_in, pt_out)
            for( n=0; n<(long long)N; n++ )
            {
                this->transform(pt_in+n*D, pt_out+n*D);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::transform(size_t* pt_in, size_t N, T* pt_out) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::transform(const size_t& xi, const size_t& yi, T& xo, T& yo) const
    {
        try
        {
            xo = xi + this->deform_field_[0](xi, yi);
            yo = yi + this->deform_field_[1](xi, yi);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::transform(const size_t& xi, const size_t& yi, T& xo, T& yo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::transform(const size_t* xi, const size_t* yi, size_t N, T* xo, T* yo) const
    {
        try
        {
            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, xi, yi, xo, yo)
            for( n=0; n<(long long)N; n++ )
            {
                this->transform(xi[n], yi[n], xo[n], yo[n]);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::transform(size_t* xi, size_t* yi, size_t N, T* xo, T* yo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::transform(const size_t& xi, const size_t& yi, const size_t& zi, T& xo, T& yo, T& zo) const
    {
        try
        {
            xo = xi + this->deform_field_[0](xi, yi, zi);
            yo = yi + this->deform_field_[1](xi, yi, zi);
            zo = zi + this->deform_field_[2](xi, yi, zi);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::transform(const size_t& xi, const size_t& yi, const size_t& zi, T& xo, T& yo, T& zo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::transform(const size_t* xi, const size_t* yi, const size_t* zi, size_t N, T* xo, T* yo, T* zo) const
    {
        try
        {
            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, xi, yi, zi, xo, yo, zo)
            for( n=0; n<(long long)N; n++ )
            {
                this->transform(xi[n], yi[n], zi[n], xo[n], yo[n], zo[n]);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::transform(size_t* xi, size_t* yi, size_t* zi, size_t N, T* xo, T* yo, T* zo) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::jacobianPosition(const input_point_type& pos, jacobian_position_type& jac)
    {
        try
        {
            jac.createMatrix(D, D);

            T delta = 0.5;
            T deltaReciprocal = T(1.0)/(T(2.0)*delta);

            std::vector<coord_type> pixelSize(D);

            coord_type pos_positive_vec[D];
            coord_type pos_negative_vec[D];

            deform_field_[0].get_pixel_size(pixelSize);

            size_t din, dout;
            for ( dout=0; dout<D; dout++ )
            {
                for ( din=0; din<D; din++ )
                {
                    input_point_type pos_positive(pos);
                    input_point_type pos_negative(pos);

                    pos_positive[din] += delta;
                    pos_negative[din] -= delta;

                    for (size_t dd = 0; dd < D; dd++)
                    {
                        pos_positive_vec[dd] = (coord_type)pos_positive[dd];
                        pos_negative_vec[dd] = (coord_type)pos_negative[dd];
                    }

                    T v_positive = (*interp_default_[dout])(pos_positive_vec);
                    T v_negative = (*interp_default_[dout])(pos_negative_vec);

                    jac(dout, din) = (v_positive-v_negative)*deltaReciprocal;

                    if ( dout != din )
                    {
                        // scaled for non-isotropic pixel sizes
                        jac(dout, din) *= ( pixelSize[dout]/pixelSize[din] );
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::jacobianPosition(const input_point_type& pos, jacobian_position_type& jac) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::jacobianPosition(hoNDArray<T>& jac, unsigned int borderWidth)
    {
        DeformationFieldType* deform_field[D];

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            deform_field[ii] = &deform_field_[ii];
        }

        return this->jacobianPosition(jac, deform_field, borderWidth);
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::jacobianPosition(hoNDArray<T>& jac, DeformationFieldType* deform_field[D], unsigned int borderWidth)
    {
        try
        {
            std::vector<size_t> dim;
            deform_field[0]->get_dimensions(dim);

            std::vector<size_t> dimJac(D+2, D);
            memcpy(&dimJac[0]+2, &dim[0], sizeof(size_t)*D);

            jac.create(dimJac);
            Gadgetron::clear(&jac);

            std::vector<size_t> offset(D);
            deform_field[0]->get_offset_factor(offset);

            std::vector<coord_type> pixelSize(D);
            deform_field[0]->get_pixel_size(pixelSize);

            T delta = 1.0;
            T deltaReciprocal = T(1.0)/(T(2.0)*delta);

            size_t N = deform_field[0]->get_number_of_elements();

            long long n;

            #pragma omp parallel private(n) shared(N, jac, dim, offset, pixelSize, borderWidth, deltaReciprocal, deform_field)
            {

                std::vector<size_t> ind(D);

                hoNDArray<T> jacCurr(D, D);

                #pragma omp for
                for ( n=0; n<(long long)N; n++ )
                {
                    ind = deform_field[0]->calculate_index( n );

                    bool inRange = true;

                    size_t din, dout;

                    for ( dout=0; dout<D; dout++ )
                    {
                        if ( ind[dout]<borderWidth || ind[dout]>=dim[dout]-borderWidth )
                        {
                            inRange = false;
                            break;
                        }
                    }

                    if ( inRange )
                    {
                        for ( dout=0; dout<D; dout++ )
                        {
                            for ( din=0; din<D; din++ )
                            {
                                size_t offset_positive = n + offset[din];
                                size_t offset_negative = n - offset[din];

                                T v_positive = (*deform_field[dout])(offset_positive);
                                T v_negative = (*deform_field[dout])(offset_negative);

                                jacCurr(dout, din) = (v_positive-v_negative)*deltaReciprocal;

                                if ( dout != din )
                                {
                                    // scaled for non-isotropic pixel sizes
                                    jacCurr(dout, din) *= ( pixelSize[dout]/pixelSize[din] );
                                }
                            }
                        }

                        memcpy(jac.begin()+n*D*D, jacCurr.begin(), sizeof(T)*D*D);
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::jacobianPosition(hoNDArray<T>& jac, DeformationFieldType* deform_field[D], unsigned int borderWidth) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline bool hoImageRegDeformationField<ValueType, D>::
    analyzeJacobianAndDeformation(const hoNDArray<T>& jac, T& meanDeform, T& maxDeform, T& meanLogJac, T& maxLogJac, unsigned int borderWidth)
    {
        DeformationFieldType* deform_field[D];

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            deform_field[ii] = &deform_field_[ii];
        }

        return this->analyzeJacobianAndDeformation(jac, deform_field, meanDeform, maxDeform, meanLogJac, maxLogJac, borderWidth);
    }

    template <typename ValueType, unsigned int D>
    bool hoImageRegDeformationField<ValueType, D>::
    analyzeJacobianAndDeformation(const hoNDArray<T>& jac, DeformationFieldType* deform_field[D], T& meanDeform, T& maxDeform, T& meanLogJac, T& maxLogJac, unsigned int borderWidth)
    {
        try
        {
            std::vector<size_t> dim;
            deform_field[0]->get_dimensions(dim);

            std::vector<coord_type> pixelSize(D);
            deform_field[0]->get_pixel_size(pixelSize);

            size_t N = deform_field[0]->get_number_of_elements();

            meanDeform = 0;
            maxDeform = -1;
            meanLogJac = 0;
            maxLogJac = -1;

            hoNDArray<T> deformNorm(dim);
            Gadgetron::clear(deformNorm);

            hoNDArray<T> logJac(dim);
            Gadgetron::clear(logJac);

            long long n;
            #pragma omp parallel private(n) shared(N, borderWidth, jac, deformNorm, logJac, dim, pixelSize, deform_field)
            {
                std::vector<size_t> ind(D);
                hoMatrix<T> jacCurr(D, D);
                unsigned int ii;

                #pragma omp for
                for ( n=0; n<(long long)N; n++ )
                {
                    ind = deform_field[0]->calculate_index( n );

                    bool inRange = true;

                    size_t dout;

                    for ( dout=0; dout<D; dout++ )
                    {
                        if ( ind[dout]<borderWidth || ind[dout]>=dim[dout]-borderWidth )
                        {
                            inRange = false;
                            break;
                        }
                    }

                    if ( inRange )
                    {
                        memcpy(jacCurr.begin(), jac.begin()+n*D*D, sizeof(T)*D*D);

                        T deformMag(0), v, det;

                        for ( ii=0; ii<D; ii++ )
                        {
                            jacCurr(ii, ii) += 1.0;

                            v = (*deform_field[ii])(n)*pixelSize[ii];
                            deformMag += v*v;
                        }

                        deformNorm(n) = std::sqrt(deformMag);

                        if ( D == 2 )
                        {
                            det = jacCurr(0, 0)*jacCurr(1, 1) - jacCurr(0, 1)*jacCurr(1, 0);
                        }
                        else if ( D == 3 )
                        {
                            det = jacCurr(0, 0)*jacCurr(1, 1)*jacCurr(2, 2)
                                + jacCurr(0, 1)*jacCurr(1, 2)*jacCurr(2, 0)
                                + jacCurr(0, 2)*jacCurr(2, 1)*jacCurr(1, 0)
                                - jacCurr(0, 2)*jacCurr(1, 1)*jacCurr(2, 0)
                                - jacCurr(0, 1)*jacCurr(1, 0)*jacCurr(2, 2)
                                - jacCurr(0, 0)*jacCurr(2, 1)*jacCurr(1, 2);
                        }

                        // if ( std::abs(det) < FLT_EPSILON ) det = FLT_EPSILON;
                        if ( det < FLT_EPSILON ) det = FLT_EPSILON;
                        logJac(n) = std::log(det);
                        if( std::isnan(logJac(n)) ) logJac(n) = 0;
                    }
                }
            }

            size_t ind;
            Gadgetron::maxAbsolute(deformNorm, maxDeform, ind);
            Gadgetron::maxAbsolute(logJac, maxLogJac, ind);

            double totalDeform = 0;
            for ( n=0; n<(long long)N; n++ )
            {
                totalDeform += deformNorm(n);
            }

            // Gadgetron::norm1(deformNorm, meanDeform);
            meanDeform = (T)(totalDeform/N);

            double totalLogJac = 0;
            for ( n=0; n<(long long)N; n++ )
            {
                totalLogJac += std::abs(logJac(n));
            }

            // Gadgetron::norm1(logJac, meanLogJac);
            meanLogJac = (T)(totalLogJac/N);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in analyzeJacobianAndDeformation(const hoNDArray<T>& jac, DeformationFieldType* deform_field[D], T& meanDeform, T& maxDeform, T& meanLogJac, T& maxLogJac, unsigned int borderWidth) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    inline ValueType& hoImageRegDeformationField<ValueType, D>::operator()( size_t idx[D], size_t outDim )
    {
        GADGET_DEBUG_CHECK_THROW(outDim<=D);
        return this->deform_field_[outDim](idx);
    }

    template <typename ValueType, unsigned int D>
    inline const ValueType& hoImageRegDeformationField<ValueType, D>::operator()( size_t idx[D], size_t outDim ) const
    {
        GADGET_DEBUG_CHECK_THROW(outDim<=D);
        return this->deform_field_[outDim](idx);
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::get(size_t idx[D], T deform[D])
    {
        size_t offset = this->deform_field_[0].calculate_offset(idx);

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            deform[ii] = this->deform_field_[ii](offset);
        }
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::get(size_t x, size_t y, T& dx, T& dy)
    {
        size_t offset = this->deform_field_[0].calculate_offset(x, y);
        dx = this->deform_field_[0](offset);
        dy = this->deform_field_[1](offset);
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::get(size_t x, size_t y, size_t z, T& dx, T& dy, T& dz)
    {
        size_t offset = this->deform_field_[0].calculate_offset(x, y, z);
        dx = this->deform_field_[0](offset);
        dy = this->deform_field_[1](offset);
        dz = this->deform_field_[2](offset);
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::set(size_t idx[D], T deform[D])
    {
        size_t offset = this->deform_field_[0].calculate_offset(idx);

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->deform_field_[ii](offset) = deform[ii];
        }
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::set(size_t x, size_t y, T dx, T dy)
    {
        size_t offset = this->deform_field_[0].calculate_offset(x, y);
        this->deform_field_[0](offset) = dx;
        this->deform_field_[1](offset) = dy;
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::set(size_t x, size_t y, size_t z, T dx, T dy, T dz)
    {
        size_t offset = this->deform_field_[0].calculate_offset(x, y, z);
        this->deform_field_[0](offset) = dx;
        this->deform_field_[1](offset) = dy;
        this->deform_field_[2](offset) = dz;
    }

    template <typename ValueType, unsigned int D>
    inline ValueType hoImageRegDeformationField<ValueType, D>::operator()( coord_type pos[D], size_t outDim )
    {
        GADGET_DEBUG_CHECK_THROW(outDim<=D);
        return (*interp_default_[outDim])(pos);
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::get(coord_type pos[D], T deform[D])
    {
        unsigned int ii;
        for (ii=0; ii<D; ii++ )
        {
            deform[ii] = (*interp_default_[ii])(pos);
        }
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::get(coord_type px, coord_type py, T& dx, T& dy)
    {
        dx = (*interp_default_[0])(px, py);
        dy = (*interp_default_[1])(px, py);
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::get(coord_type px, coord_type py, coord_type pz, T& dx, T& dy, T& dz)
    {
        dx = (*interp_default_[0])(px, py, pz);
        dy = (*interp_default_[1])(px, py, pz);
        dz = (*interp_default_[2])(px, py, pz);
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::  getDeformationField(DeformationFieldType*& deform, size_t outDim)
    {
        GADGET_DEBUG_CHECK_THROW(outDim<=D);
        deform = &(this->deform_field_[outDim]);
    }

    template <typename ValueType, unsigned int D>
    inline void hoImageRegDeformationField<ValueType, D>::setDeformationField(const DeformationFieldType& deform, size_t outDim)
    {
        GADGET_DEBUG_CHECK_THROW(outDim<=D);
        this->deform_field_[outDim] = deform;
        this->update();
    }

    template <typename ValueType, unsigned int D>
    bool hoImageRegDeformationField<ValueType, D>::serialize(char*& buf, size_t& len) const
    {
        try
        {
            if ( buf != NULL ) delete[] buf;

            char* bufInternal[D];
            size_t lenInternal[D];

            // serialize every dimension

            size_t totalLen = 0;

            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                GADGET_CHECK_RETURN_FALSE(this->deform_field_[ii].serialize(bufInternal[ii], lenInternal[ii]));
                totalLen += lenInternal[ii];
            }

            // number of dimensions + dimension vector + pixel size + origin + axis + contents
            len = sizeof(unsigned int) + totalLen;

            buf = new char[len];
            GADGET_CHECK_RETURN_FALSE(buf!=NULL);

            unsigned int NDim=D;

            size_t offset = 0;
            memcpy(buf, &NDim, sizeof(unsigned int));
            offset += sizeof(unsigned int);

            if ( NDim > 0 )
            {
                for ( ii=0; ii<D; ii++ )
                {
                    memcpy(buf+offset, bufInternal[ii], lenInternal[ii]);
                    offset += lenInternal[ii];
                }

                for ( ii=0; ii<D; ii++ )
                {
                    delete [] bufInternal[ii];
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::serialize(...) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    bool hoImageRegDeformationField<ValueType, D>::deserialize(char* buf, size_t& len)
    {
        try
        {
            unsigned int NDim;
            memcpy(&NDim, buf, sizeof(unsigned int));
            if ( NDim != D )
            {
                GERROR_STREAM("hoImageRegDeformationField<ValueType, D>::deserialize(...) : number of image dimensions does not match ... ");
                return false;
            }

            size_t offset = sizeof(unsigned int);

            unsigned int ii;

            if ( NDim > 0 )
            {
                for ( ii=0; ii<D; ii++ )
                {
                    size_t lenInternal;
                    GADGET_CHECK_RETURN_FALSE(this->deform_field_[ii].deserialize(buf+offset, lenInternal));
                    offset += lenInternal;
                }
            }

            len = offset;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationField<ValueType, D>::deserialize(...) ... ");
            return false;
        }

        return true;
    }

    template <typename ValueType, unsigned int D>
    void hoImageRegDeformationField<ValueType, D>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron deformation field geometry transformation -------------" << endl;
        os << "Deformation field dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(T).name());
        os << "Transformation data type is : " << elemTypeName << std::endl;
    }

    template <typename ValueType, unsigned int D>
    std::string hoImageRegDeformationField<ValueType, D>::transformationName() const
    {
        return std::string("hoImageRegDeformationField");
    }
}
#endif // hoImageRegDeformationField_H_
