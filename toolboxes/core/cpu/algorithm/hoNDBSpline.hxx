/** \file       hoNDBSpline.hxx
    \brief      Implementation of N-dimensional BSpline class for gadgetron
    \author     Hui Xue
*/

#include "hoNDBSpline.h"

namespace Gadgetron
{
    template <typename T, unsigned int D, typename coord_type>
    bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const hoNDArray<T>& data, unsigned int SplineDegree, hoNDArray<T>& coeff)
    {
        size_t NDim = data.get_number_of_dimensions();

        if (!coeff.dimensions_equal(data))
        {
            coeff = data;
        }

        bool res;
        switch (NDim)
        {
        case 1:
            res = this->computeBSplineCoefficients(data.begin(), data.get_size(0), SplineDegree, coeff.begin());
            break;

        case 2:
            res = this->computeBSplineCoefficients(data.begin(), data.get_size(0), data.get_size(1), SplineDegree, coeff.begin());
            break;

        case 3:
            res = this->computeBSplineCoefficients(data.begin(), data.get_size(0), data.get_size(1), data.get_size(2), SplineDegree, coeff.begin());
            break;

        case 4:
            res = this->computeBSplineCoefficients(data.begin(), data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3), SplineDegree, coeff.begin());
            break;

        case 5:
            res = this->computeBSplineCoefficients(data.begin(), data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3), data.get_size(4), SplineDegree, coeff.begin());
            break;

        case 6:
            res = this->computeBSplineCoefficients(data.begin(), data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3), data.get_size(4), data.get_size(5), SplineDegree, coeff.begin());
            break;

        case 7:
            res = this->computeBSplineCoefficients(data.begin(), data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3), data.get_size(4), data.get_size(5), data.get_size(6), SplineDegree, coeff.begin());
            break;

        case 8:
            res = this->computeBSplineCoefficients(data.begin(), data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3), data.get_size(4), data.get_size(5), data.get_size(6), data.get_size(7), SplineDegree, coeff.begin());
            break;

        case 9:
            res = this->computeBSplineCoefficients(data.begin(), data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3), data.get_size(4), data.get_size(5), data.get_size(6), data.get_size(7), data.get_size(8), SplineDegree, coeff.begin());
            break;

        default:
            boost::shared_ptr< std::vector<size_t> > dim = data.get_dimensions();
            res = this->computeBSplineCoefficients(data.begin(), *dim, SplineDegree, coeff.begin());
        }

        return res;
    }

template <typename T, unsigned int D, typename coord_type>
    inline bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const hoNDImage<T, D>& data, unsigned int SplineDegree, hoNDArray<T>& coeff)
    {
        std::vector<size_t> dim;
        data.get_dimensions(dim);
        hoNDArray<T> dataTmp(dim, const_cast<T*>(data.begin()), false);
        return this->computeBSplineCoefficients(dataTmp, SplineDegree, coeff);
    }

    template <typename T, unsigned int D, typename coord_type>
    bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, const std::vector<size_t>& dimension, unsigned int SplineDegree, T* coeff)
    {
        try
        {
            unsigned int NbPoles;
            bspline_float_type pole[4];
            this->Pole(pole, SplineDegree, NbPoles);

            GADGET_CHECK_RETURN_FALSE(D == dimension.size());

            hoNDArray<T> coeffBuf(const_cast<std::vector<size_t>&>(dimension), coeff, false);
            memcpy(coeff, data, sizeof(T) * coeffBuf.get_number_of_elements());

            unsigned int d;
            for (d = 0; d < D; d++)
            {
                long long ii;

                size_t len = dimension[d];
                size_t num = coeffBuf.get_number_of_elements() / len;

                size_t i;
                std::vector<size_t> dimUsed(D - 1);
                for (i = 0; i < D; i++)
                {
                    if (i < d)
                    {
                        dimUsed[i] = dimension[i];
                    }
                    else if (i > d)
                    {
                        dimUsed[i - 1] = dimension[i];
                    }
                }

                std::vector<size_t> offsetFactor(D - 1, 1);
                hoNDArray<T>::calculate_offset_factors(dimUsed, offsetFactor);

                //for( i = 0; i<D-1; i++ )
                //{
                //    size_t k = 1;
                //    for( j = 0; j < i; j++ )
                //    {
                //        k *= dimUsed[j];
                //    }
                //    offsetFactor[i] = k;
                //}

#pragma omp parallel private(ii) shared(coeff, coeffBuf, pole, NbPoles, dimension, num, len, offsetFactor, d)
                {
                    T* buf = new T[len];

                    std::vector<size_t> ind(D, 0);
                    std::vector<size_t> indUsed(D - 1, 0);

#pragma omp for 
                    for (ii = 0; ii < num; ii++)
                    {
                        if (d == 0)
                        {
                            memcpy(buf, coeff + ii * len, sizeof(T) * len);

                            this->ConvertToInterpolationCoefficients(buf, len, pole, NbPoles, DBL_EPSILON);

                            memcpy(coeff + ii * len, buf, sizeof(T) * len);
                        }
                        else
                        {
                            hoNDArray<T>::calculate_index(ii, offsetFactor, indUsed);

                            long long i;

                            //size_t offset = ii;
                            //for( i=D-2; i>=0; i-- )
                            //{
                            //    indUsed[i] = offset / offsetFactor[i];
                            //    offset %= offsetFactor[i];
                            //}

                            for (i = 0; i < D; i++)
                            {
                                if (i < d)
                                {
                                    ind[i] = indUsed[i];
                                }
                                else if (i > d)
                                {
                                    ind[i] = indUsed[i - 1];
                                }
                            }

                            for (i = 0; i < len; i++)
                            {
                                ind[d] = i;
                                buf[i] = coeffBuf(ind);
                            }

                            this->ConvertToInterpolationCoefficients(buf, len, pole, NbPoles, DBL_EPSILON);

                            for (i = 0; i < len; i++)
                            {
                                ind[d] = i;
                                coeffBuf(ind) = buf[i];
                            }
                        }
                    }

                    delete[] buf;
                }
            }

        }
        catch (...)
        {
            GERROR_STREAM("Error happened in hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, const std::vector<size_t>& dimension, unsigned int SplineDegree, T* coeff) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t len, unsigned int SplineDegree, T* coeff)
    {
        try
        {
            unsigned int NbPoles;
            bspline_float_type pole[4];
            this->Pole(pole, SplineDegree, NbPoles);

            memcpy(coeff, data, sizeof(T) * len);
            this->ConvertToInterpolationCoefficients(coeff, len, pole, NbPoles, DBL_EPSILON);
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t len, unsigned int SplineDegree, T* coeff) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D, typename coord_type>
    bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, unsigned int SplineDegree, T* coeff)
    {
        try
        {
            unsigned int NbPoles;
            bspline_float_type pole[4];
            this->Pole(pole, SplineDegree, NbPoles);

            // x
            long long y;
#pragma omp parallel default(none) private(y) shared(data, coeff, pole, NbPoles, sx, sy)
            {
                T* buf = new T[sx];

#pragma omp for 
                for (y = 0; y < sy; y++)
                {
                    memcpy(buf, data + y * sx, sizeof(T) * sx);

                    this->ConvertToInterpolationCoefficients(buf, sx, pole, NbPoles, DBL_EPSILON);

                    memcpy(coeff + y * sx, buf, sizeof(T) * sx);
                }

                delete[] buf;
            }

            // y
            long long x;
#pragma omp parallel default(none) private(x) shared(data, coeff, pole, NbPoles, sx, sy)
            {
                T* buf = new T[sy];

#pragma omp for 
                for (x = 0; x < sx; x++)
                {
                    size_t y;

                    for (y = 0; y < sy; y++)
                    {
                        buf[y] = coeff[x + y * sx];
                    }

                    this->ConvertToInterpolationCoefficients(buf, sy, pole, NbPoles, DBL_EPSILON);

                    for (y = 0; y < sy; y++)
                    {
                        coeff[x + y * sx] = buf[y];
                    }
                }

                delete[] buf;
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, unsigned int SplineDegree, T* coeff) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D, typename coord_type>
    bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, unsigned int SplineDegree, T* coeff)
    {
        try
        {
            unsigned int NbPoles;
            bspline_float_type pole[4];
            this->Pole(pole, SplineDegree, NbPoles);

            // x
            long long z;
#pragma omp parallel default(none) private(z) shared(data, coeff, pole, NbPoles, sx, sy, sz)
            {
                T* buf = new T[sx];

#pragma omp for 
                for (z = 0; z < sz; z++)
                {
                    for (size_t y = 0; y < sy; y++)
                    {
                        memcpy(buf, data + z * sx * sy + y * sx, sizeof(T) * sx);

                        this->ConvertToInterpolationCoefficients(buf, sx, pole, NbPoles, DBL_EPSILON);

                        memcpy(coeff + z * sx * sy + y * sx, buf, sizeof(T) * sx);
                    }
                }

                delete[] buf;
            }

            // y
#pragma omp parallel default(none) private(z) shared(data, coeff, pole, NbPoles, sx, sy, sz)
            {
                T* buf = new T[sy];

#pragma omp for 
                for (z = 0; z < sz; z++)
                {
                    for (size_t x = 0; x < sx; x++)
                    {
                        size_t y;

                        size_t offset = x + z * sx * sy;

                        for (y = 0; y < sy; y++)
                        {
                            buf[y] = coeff[offset + y * sx];
                        }

                        this->ConvertToInterpolationCoefficients(buf, sy, pole, NbPoles, DBL_EPSILON);

                        for (y = 0; y < sy; y++)
                        {
                            coeff[offset + y * sx] = buf[y];
                        }
                    }
                }

                delete[] buf;
            }

            // z
            long long x;
#pragma omp parallel default(none) private(x) shared(data, coeff, pole, NbPoles, sx, sy, sz)
            {
                T* buf = new T[sz];

#pragma omp for 
                for (x = 0; x < sx; x++)
                {
                    for (size_t y = 0; y < sy; y++)
                    {
                        size_t z;
                        size_t offset = x + y * sx;

                        for (z = 0; z < sz; z++)
                        {
                            buf[z] = coeff[offset + z * sx * sy];
                        }

                        this->ConvertToInterpolationCoefficients(buf, sz, pole, NbPoles, DBL_EPSILON);

                        for (z = 0; z < sz; z++)
                        {
                            coeff[offset + z * sx * sy] = buf[z];
                        }
                    }
                }

                delete[] buf;
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, unsigned int SplineDegree, T* coeff) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D, typename coord_type>
    bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, unsigned int SplineDegree, T* coeff)
    {
        try
        {
            unsigned int NbPoles;
            bspline_float_type pole[4];
            this->Pole(pole, SplineDegree, NbPoles);

            long long x, y, z, t;

            // x
#pragma omp parallel default(none) private(y, z, t) shared(data, coeff, pole, NbPoles, sx, sy, sz, st)
            {
                T* buf = new T[sx];

#pragma omp for 
                for (t = 0; t < st; t++)
                {
                    for (z = 0; z < sz; z++)
                    {
                        for (y = 0; y < sy; y++)
                        {
                            memcpy(buf, data + t * sx * sy * sz + z * sx * sy + y * sx, sizeof(T) * sx);

                            this->ConvertToInterpolationCoefficients(buf, sx, pole, NbPoles, DBL_EPSILON);

                            memcpy(coeff + t * sx * sy * sz + z * sx * sy + y * sx, buf, sizeof(T) * sx);
                        }
                    }
                }

                delete[] buf;
            }

            // y
#pragma omp parallel default(none) private(x, y, z, t) shared(data, coeff, pole, NbPoles, sx, sy, sz, st)
            {
                T* buf = new T[sy];

#pragma omp for 
                for (t = 0; t < st; t++)
                {
                    for (z = 0; z < sz; z++)
                    {
                        for (x = 0; x < sx; x++)
                        {
                            size_t offset = x + z * sx * sy + t * sx * sy * sz;

                            for (y = 0; y < sy; y++)
                            {
                                buf[y] = coeff[offset + y * sx];
                            }

                            this->ConvertToInterpolationCoefficients(buf, sy, pole, NbPoles, DBL_EPSILON);

                            for (y = 0; y < sy; y++)
                            {
                                coeff[offset + y * sx] = buf[y];
                            }
                        }
                    }
                }

                delete[] buf;
            }

            // z
#pragma omp parallel default(none) private(x, y, z, t) shared(data, coeff, pole, NbPoles, sx, sy, sz, st)
            {
                T* buf = new T[sz];

#pragma omp for 
                for (t = 0; t < st; t++)
                {
                    for (x = 0; x < sx; x++)
                    {
                        for (y = 0; y < sy; y++)
                        {
                            size_t offset = x + y * sx + t * sx * sy * sz;

                            for (z = 0; z < sz; z++)
                            {
                                buf[z] = coeff[offset + z * sx * sy];
                            }

                            this->ConvertToInterpolationCoefficients(buf, sz, pole, NbPoles, DBL_EPSILON);

                            for (z = 0; z < sz; z++)
                            {
                                coeff[offset + z * sx * sy] = buf[z];
                            }
                        }
                    }
                }

                delete[] buf;
            }

            // t
#pragma omp parallel default(none) private(x, y, z, t) shared(data, coeff, pole, NbPoles, sx, sy, sz, st)
            {
                T* buf = new T[st];

#pragma omp for 
                for (x = 0; x < sx; x++)
                {
                    for (y = 0; y < sy; y++)
                    {
                        for (z = 0; z < sz; z++)
                        {
                            size_t offset = x + y * sx + z * sx * sy;

                            for (t = 0; t < st; t++)
                            {
                                buf[t] = coeff[offset + t * sx * sy * sz];
                            }

                            this->ConvertToInterpolationCoefficients(buf, st, pole, NbPoles, DBL_EPSILON);

                            for (t = 0; t < st; t++)
                            {
                                coeff[offset + t * sx * sy * sz] = buf[t];
                            }
                        }
                    }
                }

                delete[] buf;
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, unsigned int SplineDegree, T* coeff) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D, typename coord_type>
    bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, unsigned int SplineDegree, T* coeff)
    {
        try
        {
            unsigned int NbPoles;
            bspline_float_type pole[4];
            this->Pole(pole, SplineDegree, NbPoles);

            long long x, y, z, t, p;

            // x
#pragma omp parallel default(none) private(y, z, t, p) shared(data, coeff, pole, NbPoles, sx, sy, sz, st, sp)
            {
                T* buf = new T[sx];

#pragma omp for 
                for (p = 0; p < sp; p++)
                {
                    for (t = 0; t < st; t++)
                    {
                        for (z = 0; z < sz; z++)
                        {
                            for (y = 0; y < sy; y++)
                            {
                                memcpy(buf, data + p * sx * sy * sz * st + t * sx * sy * sz + z * sx * sy + y * sx, sizeof(T) * sx);

                                this->ConvertToInterpolationCoefficients(buf, sx, pole, NbPoles, DBL_EPSILON);

                                memcpy(coeff + p * sx * sy * sz * st + t * sx * sy * sz + z * sx * sy + y * sx, buf, sizeof(T) * sx);
                            }
                        }
                    }
                }

                delete[] buf;
            }

            // y
#pragma omp parallel default(none) private(x, y, z, t, p) shared(data, coeff, pole, NbPoles, sx, sy, sz, st, sp)
            {
                T* buf = new T[sy];

#pragma omp for 
                for (p = 0; p < sp; p++)
                {
                    for (t = 0; t < st; t++)
                    {
                        for (z = 0; z < sz; z++)
                        {
                            for (x = 0; x < sx; x++)
                            {
                                size_t offset = x + z * sx * sy + t * sx * sy * sz + p * sx * sy * sz * st;

                                for (y = 0; y < sy; y++)
                                {
                                    buf[y] = coeff[offset + y * sx];
                                }

                                this->ConvertToInterpolationCoefficients(buf, sy, pole, NbPoles, DBL_EPSILON);

                                for (y = 0; y < sy; y++)
                                {
                                    coeff[offset + y * sx] = buf[y];
                                }
                            }
                        }
                    }
                }

                delete[] buf;
            }

            // z
#pragma omp parallel default(none) private(x, y, z, t, p) shared(data, coeff, pole, NbPoles, sx, sy, sz, st, sp)
            {
                T* buf = new T[sz];

#pragma omp for 
                for (p = 0; p < sp; p++)
                {
                    for (t = 0; t < st; t++)
                    {
                        for (x = 0; x < sx; x++)
                        {
                            for (y = 0; y < sy; y++)
                            {
                                size_t offset = x + y * sx + t * sx * sy * sz + p * sx * sy * sz * st;

                                for (z = 0; z < sz; z++)
                                {
                                    buf[z] = coeff[offset + z * sx * sy];
                                }

                                this->ConvertToInterpolationCoefficients(buf, sz, pole, NbPoles, DBL_EPSILON);

                                for (z = 0; z < sz; z++)
                                {
                                    coeff[offset + z * sx * sy] = buf[z];
                                }
                            }
                        }
                    }
                }

                delete[] buf;
            }

            // t
#pragma omp parallel default(none) private(x, y, z, t, p) shared(data, coeff, pole, NbPoles, sx, sy, sz, st, sp)
            {
                T* buf = new T[st];

#pragma omp for 
                for (p = 0; p < sp; p++)
                {
                    for (x = 0; x < sx; x++)
                    {
                        for (y = 0; y < sy; y++)
                        {
                            for (z = 0; z < sz; z++)
                            {
                                size_t offset = x + y * sx + z * sx * sy + p * sx * sy * sz * st;

                                for (t = 0; t < st; t++)
                                {
                                    buf[t] = coeff[offset + t * sx * sy * sz];
                                }

                                this->ConvertToInterpolationCoefficients(buf, st, pole, NbPoles, DBL_EPSILON);

                                for (t = 0; t < st; t++)
                                {
                                    coeff[offset + t * sx * sy * sz] = buf[t];
                                }
                            }
                        }
                    }
                }

                delete[] buf;
            }

            // p
#pragma omp parallel default(none) private(x, y, z, t, p) shared(data, coeff, pole, NbPoles, sx, sy, sz, st, sp)
            {
                T* buf = new T[sp];

#pragma omp for 
                for (x = 0; x < sx; x++)
                {
                    for (y = 0; y < sy; y++)
                    {
                        for (z = 0; z < sz; z++)
                        {
                            for (t = 0; t < st; t++)
                            {
                                size_t offset = x + y * sx + z * sx * sy + t * sx * sy * sz;

                                for (p = 0; p < sp; p++)
                                {
                                    buf[t] = coeff[offset + p * sx * sy * sz * st];
                                }

                                this->ConvertToInterpolationCoefficients(buf, sp, pole, NbPoles, DBL_EPSILON);

                                for (p = 0; p < sp; p++)
                                {
                                    coeff[offset + p * sx * sy * sz * st] = buf[t];
                                }
                            }
                        }
                    }
                }

                delete[] buf;
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, unsigned int SplineDegree, T* coeff) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, unsigned int SplineDegree, T* coeff)
    {
        std::vector<size_t> dim(6);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        dim[5] = sq;

        return this->computeBSplineCoefficients(data, dim, SplineDegree, coeff);
    }

    template <typename T, unsigned int D, typename coord_type>
    inline bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, unsigned int SplineDegree, T* coeff)
    {
        std::vector<size_t> dim(7);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        dim[5] = sq;
        dim[6] = sr;

        return this->computeBSplineCoefficients(data, dim, SplineDegree, coeff);
    }

    template <typename T, unsigned int D, typename coord_type>
    inline bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, unsigned int SplineDegree, T* coeff)
    {
        std::vector<size_t> dim(8);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        dim[5] = sq;
        dim[6] = sr;
        dim[7] = ss;

        return this->computeBSplineCoefficients(data, dim, SplineDegree, coeff);
    }

    template <typename T, unsigned int D, typename coord_type>
    inline bool hoNDBSpline<T, D, coord_type>::computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, size_t su, unsigned int SplineDegree, T* coeff)
    {
        std::vector<size_t> dim(9);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        dim[5] = sq;
        dim[6] = sr;
        dim[7] = ss;
        dim[8] = su;

        return this->computeBSplineCoefficients(data, dim, SplineDegree, coeff);
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, const std::vector<size_t>& dimension, unsigned int SplineDegree,
        const std::vector<unsigned int>& derivative, const coord_type* pos)
    {
        if (D != dimension.size())
        {
            GERROR_STREAM("D!=dimension.get_number_of_dimensions()");
            return T(0);
        }

        bspline_float_type weight[D][10];
        long long index[D][10];

        unsigned int ii, jj;
        for (ii = 0; ii < D; ii++)
        {
            computeBSplineInterpolationLocationsAndWeights(dimension[ii], SplineDegree, derivative[ii], pos[ii], weight[ii], index[ii]);
        }

        std::vector<size_t> splineDimension(D, SplineDegree);
        std::vector<size_t> splineInd(D, 0);
        std::vector<size_t> coeffInd(D, 0);

        std::vector<size_t> offsetFactors(D, 0);
        hoNDArray<T>::calculate_offset_factors(splineDimension, offsetFactors);

        std::vector<size_t> coeffOffsetFactors(D, 0);
        hoNDArray<T>::calculate_offset_factors(dimension, coeffOffsetFactors);

        unsigned int num = (unsigned int)std::pow((double)SplineDegree, (double)D);

        T res = 0;

        for (ii = 0; ii < num; ii++)
        {
            hoNDArray<T>::calculate_index(ii, offsetFactors, splineInd);

            for (jj = 0; jj < D; jj++)
            {
                coeffInd[jj] = index[jj][splineInd[jj]];
            }

            size_t offset = hoNDArray<T>::calculate_offset(coeffInd, coeffOffsetFactors);

            T v = coeff[offset];

            for (jj = 0; jj < D; jj++)
            {
                v *= weight[jj][splineInd[jj]];
            }

            res += v;
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, const std::vector<size_t>& dimension, unsigned int SplineDegree,
        const std::vector<unsigned int>& derivative, const std::vector<coord_type>& pos)
    {
        return this->evaluateBSpline(coeff, dimension, SplineDegree, derivative, &pos[0]);
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t len, unsigned int SplineDegree,
        unsigned int dx,
        coord_type x)
    {
        bspline_float_type xWeight[10];
        long long xIndex[10];

        computeBSplineInterpolationLocationsAndWeights(len, SplineDegree, dx, x, xWeight, xIndex);

        T res = 0;
        unsigned int ix;
        for (ix = 0; ix < SplineDegree; ix++)
        {
            res += coeff[xIndex[ix]] * xWeight[ix];
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, unsigned int SplineDegree,
        unsigned int dx, unsigned int dy,
        coord_type x, coord_type y)
    {
        bspline_float_type xWeight[10];
        long long xIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sx, SplineDegree, dx, x, xWeight, xIndex);

        bspline_float_type yWeight[10];
        long long yIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sy, SplineDegree, dy, y, yWeight, yIndex);

        T res = 0;

        unsigned int ix, iy;
        for (iy = 0; iy <= SplineDegree; iy++)
        {
            for (ix = 0; ix <= SplineDegree; ix++)
            {
                res += coeff[xIndex[ix] + sx * yIndex[iy]] * xWeight[ix] * yWeight[iy];
            }
        }

        return res;
    }


    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, unsigned int SplineDegree,
        unsigned int dx, unsigned int dy, unsigned int dz,
        coord_type x, coord_type y, coord_type z)
    {
        bspline_float_type xWeight[10];
        long long xIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sx, SplineDegree, dx, x, xWeight, xIndex);

        bspline_float_type yWeight[10];
        long long yIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sy, SplineDegree, dy, y, yWeight, yIndex);

        bspline_float_type zWeight[10];
        long long zIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sz, SplineDegree, dz, z, zWeight, zIndex);

        T res = 0;

        unsigned int ix, iy, iz;
        for (iz = 0; iz <= SplineDegree; iz++)
        {
            for (iy = 0; iy <= SplineDegree; iy++)
            {
                long long offset = yIndex[iy] * sx + zIndex[iz] * sx * sy;

                for (ix = 0; ix <= SplineDegree; ix++)
                {
                    res += coeff[xIndex[ix] + offset]
                        * xWeight[ix] * yWeight[iy] * zWeight[iz];
                }
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, unsigned int SplineDegree,
        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt,
        coord_type x, coord_type y, coord_type z, coord_type t)
    {
        bspline_float_type xWeight[10];
        long long xIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sx, SplineDegree, dx, x, xWeight, xIndex);

        bspline_float_type yWeight[10];
        long long yIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sy, SplineDegree, dy, y, yWeight, yIndex);

        bspline_float_type zWeight[10];
        long long zIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sz, SplineDegree, dz, z, zWeight, zIndex);

        bspline_float_type tWeight[10];
        long long tIndex[10];
        computeBSplineInterpolationLocationsAndWeights(st, SplineDegree, dt, t, tWeight, tIndex);

        T res = 0;

        unsigned int ix, iy, iz, it;
        for (it = 0; it <= SplineDegree; it++)
        {
            for (iz = 0; iz <= SplineDegree; iz++)
            {
                for (iy = 0; iy <= SplineDegree; iy++)
                {
                    long long offset = yIndex[iy] * sx + zIndex[iz] * sx * sy + tIndex[it] * sx * sy * sz;

                    for (ix = 0; ix <= SplineDegree; ix++)
                    {
                        res += coeff[xIndex[ix] + offset]
                            * xWeight[ix] * yWeight[iy] * zWeight[iz] * tWeight[it];
                    }
                }
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, unsigned int SplineDegree,
        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp,
        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p)
    {
        bspline_float_type xWeight[10];
        long long xIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sx, SplineDegree, dx, x, xWeight, xIndex);

        bspline_float_type yWeight[10];
        long long yIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sy, SplineDegree, dy, y, yWeight, yIndex);

        bspline_float_type zWeight[10];
        long long zIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sz, SplineDegree, dz, z, zWeight, zIndex);

        bspline_float_type tWeight[10];
        long long tIndex[10];
        computeBSplineInterpolationLocationsAndWeights(st, SplineDegree, dt, t, tWeight, tIndex);

        bspline_float_type pWeight[10];
        long long pIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sp, SplineDegree, dp, p, pWeight, pIndex);

        T res = 0;

        unsigned int ix, iy, iz, it, ip;

        for (ip = 0; ip <= SplineDegree; ip++)
        {
            for (it = 0; it <= SplineDegree; it++)
            {
                for (iz = 0; iz <= SplineDegree; iz++)
                {
                    for (iy = 0; iy <= SplineDegree; iy++)
                    {
                        long long offset = yIndex[iy] * sx + zIndex[iz] * sx * sy + tIndex[it] * sx * sy * sz + pIndex[ip] * sx * sy * sz * st;

                        for (ix = 0; ix <= SplineDegree; ix++)
                        {
                            res += coeff[xIndex[ix] + offset]
                                * xWeight[ix] * yWeight[iy] * zWeight[iz] * tWeight[it] * pWeight[ip];
                        }
                    }
                }
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, unsigned int SplineDegree,
        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp, unsigned int dq,
        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p, coord_type q)
    {
        bspline_float_type xWeight[10];
        long long xIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sx, SplineDegree, dx, x, xWeight, xIndex);

        bspline_float_type yWeight[10];
        long long yIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sy, SplineDegree, dy, y, yWeight, yIndex);

        bspline_float_type zWeight[10];
        long long zIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sz, SplineDegree, dz, z, zWeight, zIndex);

        bspline_float_type tWeight[10];
        long long tIndex[10];
        computeBSplineInterpolationLocationsAndWeights(st, SplineDegree, dt, t, tWeight, tIndex);

        bspline_float_type pWeight[10];
        long long pIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sp, SplineDegree, dp, p, pWeight, pIndex);

        bspline_float_type qWeight[10];
        long long qIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sq, SplineDegree, dq, q, qWeight, qIndex);

        T res = 0;

        unsigned int ix, iy, iz, it, ip, iq;

        for (iq = 0; iq <= SplineDegree; iq++)
        {
            for (ip = 0; ip <= SplineDegree; ip++)
            {
                for (it = 0; it <= SplineDegree; it++)
                {
                    for (iz = 0; iz <= SplineDegree; iz++)
                    {
                        for (iy = 0; iy <= SplineDegree; iy++)
                        {
                            long long offset = yIndex[iy] * sx
                                + zIndex[iz] * sx * sy
                                + tIndex[it] * sx * sy * sz
                                + pIndex[ip] * sx * sy * sz * st
                                + qIndex[iq] * sx * sy * sz * st * sp;

                            for (ix = 0; ix <= SplineDegree; ix++)
                            {
                                res += coeff[xIndex[ix] + offset]
                                    * xWeight[ix] * yWeight[iy] * zWeight[iz] * tWeight[it] * pWeight[ip] * qWeight[iq];
                            }
                        }
                    }
                }
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, unsigned int SplineDegree,
        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp, unsigned int dq, unsigned int dr,
        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p, coord_type q, coord_type r)
    {
        bspline_float_type xWeight[10];
        long long xIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sx, SplineDegree, dx, x, xWeight, xIndex);

        bspline_float_type yWeight[10];
        long long yIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sy, SplineDegree, dy, y, yWeight, yIndex);

        bspline_float_type zWeight[10];
        long long zIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sz, SplineDegree, dz, z, zWeight, zIndex);

        bspline_float_type tWeight[10];
        long long tIndex[10];
        computeBSplineInterpolationLocationsAndWeights(st, SplineDegree, dt, t, tWeight, tIndex);

        bspline_float_type pWeight[10];
        long long pIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sp, SplineDegree, dp, p, pWeight, pIndex);

        bspline_float_type qWeight[10];
        long long qIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sq, SplineDegree, dq, q, qWeight, qIndex);

        bspline_float_type rWeight[10];
        long long rIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sr, SplineDegree, dr, r, rWeight, rIndex);

        T res = 0;

        unsigned int ix, iy, iz, it, ip, iq, ir;

        for (ir = 0; ir <= SplineDegree; ir++)
        {
            for (iq = 0; iq <= SplineDegree; iq++)
            {
                for (ip = 0; ip <= SplineDegree; ip++)
                {
                    for (it = 0; it <= SplineDegree; it++)
                    {
                        for (iz = 0; iz <= SplineDegree; iz++)
                        {
                            for (iy = 0; iy <= SplineDegree; iy++)
                            {
                                long long offset = yIndex[iy] * sx
                                    + zIndex[iz] * sx * sy
                                    + tIndex[it] * sx * sy * sz
                                    + pIndex[ip] * sx * sy * sz * st
                                    + qIndex[iq] * sx * sy * sz * st * sp
                                    + rIndex[ir] * sx * sy * sz * st * sp * sq;

                                for (ix = 0; ix <= SplineDegree; ix++)
                                {
                                    res += coeff[xIndex[ix] + offset]
                                        * xWeight[ix] * yWeight[iy] * zWeight[iz] * tWeight[it] * pWeight[ip] * qWeight[iq] * rWeight[ir];
                                }
                            }
                        }
                    }
                }
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, unsigned int SplineDegree,
        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp, unsigned int dq, unsigned int dr, unsigned int ds,
        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p, coord_type q, coord_type r, coord_type s)
    {
        bspline_float_type xWeight[10];
        long long xIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sx, SplineDegree, dx, x, xWeight, xIndex);

        bspline_float_type yWeight[10];
        long long yIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sy, SplineDegree, dy, y, yWeight, yIndex);

        bspline_float_type zWeight[10];
        long long zIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sz, SplineDegree, dz, z, zWeight, zIndex);

        bspline_float_type tWeight[10];
        long long tIndex[10];
        computeBSplineInterpolationLocationsAndWeights(st, SplineDegree, dt, t, tWeight, tIndex);

        bspline_float_type pWeight[10];
        long long pIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sp, SplineDegree, dp, p, pWeight, pIndex);

        bspline_float_type qWeight[10];
        long long qIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sq, SplineDegree, dq, q, qWeight, qIndex);

        bspline_float_type rWeight[10];
        long long rIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sr, SplineDegree, dr, r, rWeight, rIndex);

        bspline_float_type sWeight[10];
        long long sIndex[10];
        computeBSplineInterpolationLocationsAndWeights(ss, SplineDegree, ds, s, sWeight, sIndex);

        T res = 0;

        unsigned int ix, iy, iz, it, ip, iq, ir, is;

        for (is = 0; is <= SplineDegree; is++)
        {
            for (ir = 0; ir <= SplineDegree; ir++)
            {
                for (iq = 0; iq <= SplineDegree; iq++)
                {
                    for (ip = 0; ip <= SplineDegree; ip++)
                    {
                        for (it = 0; it <= SplineDegree; it++)
                        {
                            for (iz = 0; iz <= SplineDegree; iz++)
                            {
                                for (iy = 0; iy <= SplineDegree; iy++)
                                {
                                    long long offset = yIndex[iy] * sx
                                        + zIndex[iz] * sx * sy
                                        + tIndex[it] * sx * sy * sz
                                        + pIndex[ip] * sx * sy * sz * st
                                        + qIndex[iq] * sx * sy * sz * st * sp
                                        + rIndex[ir] * sx * sy * sz * st * sp * sq
                                        + sIndex[ir] * sx * sy * sz * st * sp * sq * sr;

                                    for (ix = 0; ix <= SplineDegree; ix++)
                                    {
                                        res += coeff[xIndex[ix] + offset]
                                            * xWeight[ix] * yWeight[iy] * zWeight[iz] * tWeight[it] * pWeight[ip] * qWeight[iq] * rWeight[ir] * sWeight[is];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, size_t su, unsigned int SplineDegree,
        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp, unsigned int dq, unsigned int dr, unsigned int ds, unsigned int du,
        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p, coord_type q, coord_type r, coord_type s, coord_type u)
    {
        bspline_float_type xWeight[10];
        long long xIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sx, SplineDegree, dx, x, xWeight, xIndex);

        bspline_float_type yWeight[10];
        long long yIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sy, SplineDegree, dy, y, yWeight, yIndex);

        bspline_float_type zWeight[10];
        long long zIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sz, SplineDegree, dz, z, zWeight, zIndex);

        bspline_float_type tWeight[10];
        long long tIndex[10];
        computeBSplineInterpolationLocationsAndWeights(st, SplineDegree, dt, t, tWeight, tIndex);

        bspline_float_type pWeight[10];
        long long pIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sp, SplineDegree, dp, p, pWeight, pIndex);

        bspline_float_type qWeight[10];
        long long qIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sq, SplineDegree, dq, q, qWeight, qIndex);

        bspline_float_type rWeight[10];
        long long rIndex[10];
        computeBSplineInterpolationLocationsAndWeights(sr, SplineDegree, dr, r, rWeight, rIndex);

        bspline_float_type sWeight[10];
        long long sIndex[10];
        computeBSplineInterpolationLocationsAndWeights(ss, SplineDegree, ds, s, sWeight, sIndex);

        bspline_float_type uWeight[10];
        long long uIndex[10];
        computeBSplineInterpolationLocationsAndWeights(su, SplineDegree, du, u, uWeight, uIndex);

        T res = 0;

        unsigned int ix, iy, iz, it, ip, iq, ir, is, iu;

        for (iu = 0; iu <= SplineDegree; iu++)
        {
            for (is = 0; is <= SplineDegree; is++)
            {
                for (ir = 0; ir <= SplineDegree; ir++)
                {
                    for (iq = 0; iq <= SplineDegree; iq++)
                    {
                        for (ip = 0; ip <= SplineDegree; ip++)
                        {
                            for (it = 0; it <= SplineDegree; it++)
                            {
                                for (iz = 0; iz <= SplineDegree; iz++)
                                {
                                    for (iy = 0; iy <= SplineDegree; iy++)
                                    {
                                        long long offset = yIndex[iy] * sx
                                            + zIndex[iz] * sx * sy
                                            + tIndex[it] * sx * sy * sz
                                            + pIndex[ip] * sx * sy * sz * st
                                            + qIndex[iq] * sx * sy * sz * st * sp
                                            + rIndex[ir] * sx * sy * sz * st * sp * sq
                                            + sIndex[ir] * sx * sy * sz * st * sp * sq * sr
                                            + uIndex[ir] * sx * sy * sz * st * sp * sq * sr * ss;

                                        for (ix = 0; ix <= SplineDegree; ix++)
                                        {
                                            res += coeff[xIndex[ix] + offset]
                                                * xWeight[ix] * yWeight[iy] * zWeight[iz] * tWeight[it] * pWeight[ip] * qWeight[iq] * rWeight[ir] * sWeight[is] * uWeight[iu];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, unsigned int SplineDegree,
        bspline_float_type* xWeight, bspline_float_type* yWeight,
        coord_type x, coord_type y)
    {
        long long xIndex[10];
        BSplineInterpolationLocation(x, SplineDegree, xIndex);
        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, xIndex, sx);

        long long yIndex[10];
        BSplineInterpolationLocation(y, SplineDegree, yIndex);
        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, yIndex, sy);

        T res = 0;

        unsigned int ix, iy;
        for (iy = 0; iy <= SplineDegree; iy++)
        {
            for (ix = 0; ix <= SplineDegree; ix++)
            {
                res += coeff[xIndex[ix] + sx * yIndex[iy]] * xWeight[ix] * yWeight[iy];
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, unsigned int SplineDegree,
        bspline_float_type* xWeight, bspline_float_type* yWeight, bspline_float_type* zWeight,
        coord_type x, coord_type y, coord_type z)
    {
        long long xIndex[10];
        BSplineInterpolationLocation(x, SplineDegree, xIndex);
        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, xIndex, sx);

        long long yIndex[10];
        BSplineInterpolationLocation(y, SplineDegree, yIndex);
        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, yIndex, sy);

        long long zIndex[10];
        BSplineInterpolationLocation(z, SplineDegree, zIndex);
        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, zIndex, sz);

        T res = 0;

        unsigned int ix, iy, iz;
        for (iz = 0; iz <= SplineDegree; iz++)
        {
            for (iy = 0; iy <= SplineDegree; iy++)
            {
                long long offset = yIndex[iy] * sx + zIndex[iz] * sx * sy;

                for (ix = 0; ix <= SplineDegree; ix++)
                {
                    res += coeff[xIndex[ix] + offset]
                        * xWeight[ix] * yWeight[iy] * zWeight[iz];
                }
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, unsigned int SplineDegree,
        bspline_float_type* xWeight, bspline_float_type* yWeight, bspline_float_type* zWeight, bspline_float_type* tWeight,
        coord_type x, coord_type y, coord_type z, coord_type t)
    {
        long long xIndex[10];
        BSplineInterpolationLocation(x, SplineDegree, xIndex);
        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, xIndex, sx);

        long long yIndex[10];
        BSplineInterpolationLocation(y, SplineDegree, yIndex);
        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, yIndex, sy);

        long long zIndex[10];
        BSplineInterpolationLocation(z, SplineDegree, zIndex);
        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, zIndex, sz);

        long long tIndex[10];
        BSplineInterpolationLocation(t, SplineDegree, tIndex);
        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, tIndex, st);

        T res = 0;

        unsigned int ix, iy, iz, it;
        for (it = 0; it <= SplineDegree; it++)
        {
            for (iz = 0; iz <= SplineDegree; iz++)
            {
                for (iy = 0; iy <= SplineDegree; iy++)
                {
                    long long offset = yIndex[iy] * sx + zIndex[iz] * sx * sy + tIndex[it] * sx * sy * sz;

                    for (ix = 0; ix <= SplineDegree; ix++)
                    {
                        res += coeff[xIndex[ix] + offset]
                            * xWeight[ix] * yWeight[iy] * zWeight[iz] * tWeight[it];
                    }
                }
            }
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, const std::vector<size_t>& dimension, unsigned int SplineDegree,
        bspline_float_type** weight, const coord_type* pos)
    {
        long long index[D][10];

        unsigned int ii, jj;
        for (ii = 0; ii < D; ii++)
        {
            BSplineInterpolationLocation(pos[ii], SplineDegree, index[ii]);
            BSplineInterpolationMirrorBoundaryCondition(SplineDegree, index[ii], dimension[ii]);
        }

        std::vector<size_t> splineDimension(D, SplineDegree);
        std::vector<size_t> splineInd(D, 0);
        std::vector<size_t> coeffInd(D, 0);

        std::vector<size_t> offsetFactors(D, 0);
        hoNDArray<T>::calculate_offset_factors(splineDimension, offsetFactors);

        std::vector<size_t> coeffOffsetFactors(D, 0);
        hoNDArray<T>::calculate_offset_factors(dimension, coeffOffsetFactors);

        unsigned int num = pow(SplineDegree, D);

        T res = 0;

        for (ii = 0; ii < num; ii++)
        {
            hoNDArray<T>::calculate_index(ii, offsetFactors, splineInd);

            for (jj = 0; jj < D; jj++)
            {
                coeffInd[jj] = index[jj][splineInd[jj]];
            }

            size_t offset = hoNDArray<T>::calculate_offset(coeffInd, coeffOffsetFactors);

            T v = coeff[offset];

            for (jj = 0; jj < D; jj++)
            {
                v *= weight[jj][splineInd[jj]];
            }

            res += v;
        }

        return res;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline T hoNDBSpline<T, D, coord_type>::evaluateBSpline(const T* coeff, const std::vector<size_t>& dimension, unsigned int SplineDegree,
        bspline_float_type** weight, const std::vector<coord_type>& pos)
    {
        return this->evaluateBSpline(coeff, dimension, SplineDegree, weight, &pos[0]);
    }

    template <typename T, unsigned int D, typename coord_type>
    bool hoNDBSpline<T, D, coord_type>::computeBSplineDerivative(const hoNDArray<T>& data, const hoNDArray<T>& coeff, unsigned int SplineDegree, const std::vector<unsigned int>& derivative, hoNDArray<T>& deriv)
    {
        try
        {
            std::vector<size_t> dimension;
            data.get_dimensions(dimension);

            if (D != data.get_number_of_dimensions())
            {
                GERROR_STREAM("computeBSplineDerivative(hoNDArray) : D!=dimension.get_number_of_dimensions() ... ");
                return T(0);
            }

            if (!deriv.dimensions_equal(data))
            {
                deriv.create(data.get_dimensions());
            }

            // only need to compute the weights once, since this is the integer point computation
            bspline_float_type weight[D][10];
            long long index[D][10];

            unsigned int ii;
            for (ii = 0; ii < D; ii++)
            {
                computeBSplineInterpolationLocationsAndWeights(dimension[ii], SplineDegree, derivative[ii], dimension[ii] / 2, weight[ii], index[ii]);
            }

            if (D == 2)
            {
                size_t sx = data.get_size(0);
                size_t sy = data.get_size(1);

                long long y;

#pragma omp parallel for default(none) private(y) shared(sx, sy, deriv, coeff, SplineDegree, weight)
                for (y = 0; y < sy; y++)
                {
                    for (size_t x = 0; x < sx; x++)
                    {
                        deriv(x, y) = evaluateBSpline(coeff.begin(), sx, sy, SplineDegree, weight[0], weight[1], x, y);
                    }
                }
            }
            else if (D == 3)
            {
                size_t sx = data.get_size(0);
                size_t sy = data.get_size(1);
                size_t sz = data.get_size(2);

                long long z;

#pragma omp parallel for default(none) private(z) shared(sx, sy, sz, deriv, coeff, SplineDegree, weight)
                for (z = 0; z < sz; z++)
                {
                    for (size_t y = 0; y < sy; y++)
                    {
                        for (size_t x = 0; x < sx; x++)
                        {
                            deriv(x, y, z) = evaluateBSpline(coeff.begin(), sx, sy, sz, SplineDegree, weight[0], weight[1], weight[2], x, y, z);
                        }
                    }
                }
            }
            else if (D == 4)
            {
                size_t sx = data.get_size(0);
                size_t sy = data.get_size(1);
                size_t sz = data.get_size(2);
                size_t st = data.get_size(3);

                long long t;

#pragma omp parallel for default(none) private(t) shared(sx, sy, sz, st, deriv, coeff, SplineDegree, weight)
                for (t = 0; t < st; t++)
                {
                    for (size_t z = 0; z < sz; z++)
                    {
                        for (size_t y = 0; y < sy; y++)
                        {
                            for (size_t x = 0; x < sx; x++)
                            {
                                deriv(x, y, z, t) = evaluateBSpline(coeff.begin(), sx, sy, sz, st, SplineDegree, weight[0], weight[1], weight[2], weight[3], x, y, z, t);
                            }
                        }
                    }
                }
            }
            else
            {
                size_t num = data.get_number_of_elements();

                long long ii;
#pragma omp parallel private(ii) shared(num, data, coeff, dimension, SplineDegree, weight)
                {
                    std::vector<size_t> ind(D);
                    std::vector<coord_type> pos(D);

#pragma omp for 
                    for (ii = 0; ii < num; ii++)
                    {
                        data.calculate_index(ii, ind);

                        for (unsigned int jj = 0; jj < D; jj++)
                        {
                            pos[jj] = ind[jj];
                        }

                        deriv(ii) = evaluateBSpline(coeff.begin(), dimension, SplineDegree, (bspline_float_type * *)weight, pos);
                    }
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in hoNDBSpline<T, D, coord_type>::computeBSplineDerivative(const hoNDArray<T>& data, const hoNDArray<T>& coeff, const std::vector<unsigned int>& derivative, hoNDArray<T>& deriv) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D, typename coord_type>
    bool hoNDBSpline<T, D, coord_type>::computeBSplineDerivative(const hoNDImage<T, D>& data, const hoNDArray<T>& coeff, unsigned int SplineDegree, const std::vector<unsigned int>& derivative, hoNDImage<T, D>& deriv)
    {
        hoNDArray<T> dataTmp(data.get_dimensions(), const_cast<T*>(data.begin()), false);

        if (!deriv.dimension_equal(&data))
        {
            deriv = data;
        }

        hoNDArray<T> derivTmp(deriv.get_dimensions(), deriv.begin(), false);

        return computeBSplineDerivative(dataTmp, coeff, SplineDegree, derivative, derivTmp);
    }

    template <typename T, unsigned int D, typename coord_type>
    bool hoNDBSpline<T, D, coord_type>::computeBSplineDerivativePoints(const hoNDArray<T>& pts, const hoNDArray<T>& coeff, unsigned int SplineDegree, const std::vector<unsigned int>& derivative, hoNDArray<T>& deriv)
    {
        try
        {
            size_t N = pts.get_size(0);
            size_t pt_D = pts.get_size(1);

            if (D != pt_D)
            {
                GERROR_STREAM("computeBSplineDerivativePoints(hoNDArray) : D!=pt_D ... ");
                return T(0);
            }

            std::vector<size_t> dimension;
            coeff.get_dimensions(dimension);

            deriv.create(N, 1);

            // only need to compute the weights once, since this is the integer point computation
            bspline_float_type weight[D][10];
            long long index[D][10];

            unsigned int ii;
            for (ii = 0; ii < D; ii++)
            {
                computeBSplineInterpolationLocationsAndWeights(dimension[ii], SplineDegree, derivative[ii], dimension[ii] / 2, weight[ii], index[ii]);
            }

            long long n;

            if (D == 2)
            {
                size_t sx = coeff.get_size(0);
                size_t sy = coeff.get_size(1);

#pragma omp parallel for default(none) private(n) shared(sx, sy, deriv, coeff, SplineDegree, weight, N, pts)
                for (n = 0; n < N; n++)
                {
                    coord_type x = pts(n, 0);
                    coord_type y = pts(n, 1);
                    deriv(n, 0) = evaluateBSpline(coeff.begin(), sx, sy, SplineDegree, weight[0], weight[1], x, y);
                }
            }
            else if (D == 3)
            {
                size_t sx = coeff.get_size(0);
                size_t sy = coeff.get_size(1);
                size_t sz = coeff.get_size(2);

#pragma omp parallel for default(none) private(n) shared(sx, sy, sz, deriv, coeff, SplineDegree, weight, N, pts)
                for (n = 0; n < N; n++)
                {
                    coord_type x = pts(n, 0);
                    coord_type y = pts(n, 1);
                    coord_type z = pts(n, 2);
                    deriv(n, 0) = evaluateBSpline(coeff.begin(), sx, sy, sz, SplineDegree, weight[0], weight[1], weight[2], x, y, z);
                }
            }
            else if (D == 4)
            {
                size_t sx = coeff.get_size(0);
                size_t sy = coeff.get_size(1);
                size_t sz = coeff.get_size(2);
                size_t st = coeff.get_size(3);

#pragma omp parallel for default(none) private(n) shared(sx, sy, sz, st, deriv, coeff, SplineDegree, weight, N, pts)
                for (n = 0; n < N; n++)
                {
                    coord_type x = pts(n, 0);
                    coord_type y = pts(n, 1);
                    coord_type z = pts(n, 2);
                    coord_type t = pts(n, 3);
                    deriv(n, 0) = evaluateBSpline(coeff.begin(), sx, sy, sz, st, SplineDegree, weight[0], weight[1], weight[2], weight[3], x, y, z, t);
                }
            }
            else
            {
                long long ii;
#pragma omp parallel private(ii) shared(coeff, dimension, SplineDegree, weight, N, pts)
                {
                    std::vector<coord_type> pos(D);

#pragma omp for 
                    for (ii = 0; ii < N; ii++)
                    {
                        for (unsigned int jj = 0; jj < D; jj++)
                        {
                            pos[jj] = pts(ii, jj);
                        }

                        deriv(ii, 0) = evaluateBSpline(coeff.begin(), dimension, SplineDegree, (bspline_float_type * *)weight, pos);
                    }
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in hoNDBSpline<T, D, coord_type>::computeBSplineDerivativePoints(const hoNDArray<T>& pts, const hoNDArray<T>& coeff, const std::vector<unsigned int>& derivative, hoNDArray<T>& deriv) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D, typename coord_type>
    void hoNDBSpline<T, D, coord_type>::print(std::ostream& os) const
    {
        using namespace std;

        os << "--------------Gagdgetron ND BSpline -------------" << endl;
        os << "Dimension is : " << D << endl;
        std::string elemTypeName = std::string(typeid(T).name());
        os << "Data type is : " << elemTypeName << std::endl;
    }

    template <typename T, unsigned int D, typename coord_type>
    void hoNDBSpline<T, D, coord_type>::ConvertToInterpolationCoefficients(T c[], size_t DataLength, bspline_float_type z[], long NbPoles, bspline_float_type Tolerance)
    { /* begin ConvertToInterpolationCoefficients */

        double Lambda = 1.0;
        long n, k;

        /* special case required by mirror boundaries */
        if (DataLength == 1L)
        {
            return;
        }

        /* compute the overall gain */
        for (k = 0L; k < NbPoles; k++)
        {
            Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
        }

        /* apply the gain */
        for (n = 0L; n < DataLength; n++)
        {
            c[n] *= Lambda;
        }

        /* loop over all poles */
        for (k = 0L; k < NbPoles; k++)
        {
            /* causal initialization */
            c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);

            /* causal recursion */
            for (n = 1L; n < DataLength; n++)
            {
                c[n] += z[k] * c[n - 1L];
            }

            /* anticausal initialization */
            c[DataLength - 1L] = InitialAntiCausalCoefficient(c, DataLength, z[k]);

            /* anticausal recursion */
            for (n = DataLength - 2L; 0 <= n; n--)
            {
                c[n] = z[k] * (c[n + 1L] - c[n]);
            }
        }
    } /* end ConvertToInterpolationCoefficients */

    template <typename T, unsigned int D, typename coord_type>
    T hoNDBSpline<T, D, coord_type>::InitialCausalCoefficient(T c[], size_t DataLength, bspline_float_type z, bspline_float_type Tolerance)
    { /* begin InitialCausalCoefficient */

        T Sum;
        bspline_float_type zn, z2n, iz;
        size_t n, Horizon;

        /* this initialization corresponds to mirror boundaries */
        Horizon = DataLength;
        if (Tolerance > 0.0)
        {
            Horizon = (size_t)std::ceil(log(Tolerance) / log(fabs(z)));
        }

        if (Horizon < DataLength)
        {
            /* accelerated loop */
            zn = z;
            Sum = c[0];
            for (n = 1; n < Horizon; n++) {
                Sum += zn * c[n];
                zn *= z;
            }
            return(Sum);
        }
        else
        {
            /* full loop */
            zn = z;
            iz = (bspline_float_type)(1.0) / z;
            z2n = pow(z, (bspline_float_type)(DataLength - 1L));
            Sum = c[0] + z2n * c[DataLength - 1L];
            z2n *= z2n * iz;
            for (n = 1L; n <= DataLength - 2L; n++)
            {
                Sum += (zn + z2n) * c[n];
                zn *= z;
                z2n *= iz;
            }
            return(Sum / (bspline_float_type)(1.0 - zn * zn));
        }
    } /* end InitialCausalCoefficient */

    template <typename T, unsigned int D, typename coord_type>
    T hoNDBSpline<T, D, coord_type>::InitialAntiCausalCoefficient(T c[], size_t DataLength, bspline_float_type z)
    { /* begin InitialAntiCausalCoefficient */

        /* this initialization corresponds to mirror boundaries */
        return((z / (z * z - (bspline_float_type)1.0)) * (z * c[DataLength - 2L] + c[DataLength - 1L]));
    } /* end InitialAntiCausalCoefficient */

    template <typename T, unsigned int D, typename coord_type>
    inline void hoNDBSpline<T, D, coord_type>::Pole(bspline_float_type * Pole, unsigned int SplineDegree, unsigned int& NbPoles)
    {
        switch (SplineDegree)
        {
        case 2:
            NbPoles = 1;
            Pole[0] = (bspline_float_type)(std::sqrt(8.0) - 3.0);
            break;

        case 3:
            NbPoles = 1;
            Pole[0] = (bspline_float_type)(std::sqrt(3.0) - 2.0);
            break;

        case 4:
            NbPoles = 2;
            Pole[0] = (bspline_float_type)(std::sqrt(664.0 - std::sqrt(438976.0)) + std::sqrt(304.0) - 19.0);
            Pole[1] = (bspline_float_type)(std::sqrt(664.0 + std::sqrt(438976.0)) - std::sqrt(304.0) - 19.0);
            break;

        case 5:
            NbPoles = 2;
            Pole[0] = (bspline_float_type)(std::sqrt(135.0 / 2.0 - std::sqrt(17745.0 / 4.0)) + std::sqrt(105.0 / 4.0)
                - 13.0 / 2.0);
            Pole[1] = (bspline_float_type)(std::sqrt(135.0 / 2.0 + std::sqrt(17745.0 / 4.0)) - std::sqrt(105.0 / 4.0)
                - 13.0 / 2.0);
            break;

        case 6:
            NbPoles = 3;
            Pole[0] = (bspline_float_type)(-0.48829458930304475513011803888378906211227916123938);
            Pole[1] = (bspline_float_type)(-0.081679271076237512597937765737059080653379610398148);
            Pole[2] = (bspline_float_type)(-0.0014141518083258177510872439765585925278641690553467);
            break;

        case 7:
            NbPoles = 3;
            Pole[0] = (bspline_float_type)(-0.53528043079643816554240378168164607183392315234269);
            Pole[1] = (bspline_float_type)(-0.12255461519232669051527226435935734360548654942730);
            Pole[2] = (bspline_float_type)(-0.0091486948096082769285930216516478534156925639545994);
            break;

        case 8:
            NbPoles = 4;
            Pole[0] = (bspline_float_type)(-0.57468690924876543053013930412874542429066157804125);
            Pole[1] = (bspline_float_type)(-0.16303526929728093524055189686073705223476814550830);
            Pole[2] = (bspline_float_type)(-0.023632294694844850023403919296361320612665920854629);
            Pole[3] = (bspline_float_type)(-0.00015382131064169091173935253018402160762964054070043);
            break;

        case 9:
            NbPoles = 4;
            Pole[0] = (bspline_float_type)(-0.60799738916862577900772082395428976943963471853991);
            Pole[1] = (bspline_float_type)(-0.20175052019315323879606468505597043468089886575747);
            Pole[2] = (bspline_float_type)(-0.043222608540481752133321142979429688265852380231497);
            Pole[3] = (bspline_float_type)(-0.0021213069031808184203048965578486234220548560988624);
            break;

        default:
            GERROR_STREAM("Only 2 - 9 order BSpline is supported ... ");
            return;
        }
    }

    template <typename T, unsigned int D, typename coord_type>
    inline typename hoNDBSpline<T, D, coord_type>::bspline_float_type hoNDBSpline<T, D, coord_type>::BSpline(bspline_float_type x, unsigned int SplineDegree)
    {
        if (x < -((bspline_float_type)SplineDegree + 1) / 2.0)
        {
            return 0.0;
        }

        // follow the notation of origin paper

        unsigned int j, t;

        bspline_float_type value = 0.0;
        for (j = 0; j <= SplineDegree + 1; j++)
        {
            if ((x - j + 0.5 * (SplineDegree + 1)) >= 0)
            {
                bspline_float_type v1 = 1.0;
                for (t = 1; t <= j; t++)
                {
                    v1 *= t;
                }

                bspline_float_type v2 = 1.0;
                for (t = 1; t <= SplineDegree + 1 - j; t++)
                {
                    v2 *= t;
                }

                value += (bspline_float_type)((std::pow(double(-1), double(j)) * (SplineDegree + 1) / (v2 * v1)) * std::pow(x - j + 0.5 * (SplineDegree + 1), double(SplineDegree)));
            }
        }

        return value;
    }

    template <typename T, unsigned int D, typename coord_type>
    inline void hoNDBSpline<T, D, coord_type>::BSplineInterpolationLocation(bspline_float_type x, unsigned int SplineDegree, long long* xIndex)
    {
        long long i, k;

        /* compute the interpolation indexes */
        if (SplineDegree & 1L)
        {
            i = (long long)std::floor(x) - (long long)SplineDegree / 2L;
            for (k = 0L; k <= SplineDegree; k++)
            {
                xIndex[k] = i++;
            }
        }
        else
        {
            i = (long long)std::floor(x + 0.5) - (long long)SplineDegree / 2L;
            for (k = 0L; k <= SplineDegree; k++)
            {
                xIndex[k] = i++;
            }
        }
    }

    template <typename T, unsigned int D, typename coord_type>
    inline void hoNDBSpline<T, D, coord_type>::BSplineInterpolationMirrorBoundaryCondition(unsigned int SplineDegree, long long* xIndex, size_t Width)
    {
        long long Width2 = 2 * Width - 2;

        unsigned int k;

        /* apply the mirror boundary conditions */
        for (k = 0; k <= SplineDegree; k++)
        {
            xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
                (-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
                : (xIndex[k] - Width2 * (xIndex[k] / Width2)));

            if (Width <= xIndex[k])
            {
                xIndex[k] = Width2 - xIndex[k];
            }
        }
    }

    template <typename T, unsigned int D, typename coord_type>
    inline void hoNDBSpline<T, D, coord_type>::BSplineDiscrete(bspline_float_type x, unsigned int SplineDegree, bspline_float_type * xWeight, long long* xIndex)
    {
        bspline_float_type w, w2, w4, t, t0, t1;

        ///* compute the interpolation indexes */
        //if (SplineDegree & 1L)
        //{
        //    i = (long)std::floor(x) - SplineDegree / 2L;
        //    for (k = 0L; k <= SplineDegree; k++)
        //    {
        //        xIndex[k] = i++;
        //    }
        //}
        //else
        //{
        //    i = (long)std::floor(x + 0.5) - SplineDegree / 2L;
        //    for (k = 0L; k <= SplineDegree; k++)
        //    {
        //        xIndex[k] = i++;
        //    }
        //}

        // BSplineInterpolationLocation(x, SplineDegree, xIndex);

        /* compute the interpolation weights */
        switch (SplineDegree)
        {
        case 2L:
            /* x */
            w = x - (bspline_float_type)xIndex[1];
            xWeight[1] = 3.0 / 4.0 - w * w;
            xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
            xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
            break;
        case 3L:
            /* x */
            w = x - (bspline_float_type)xIndex[1];
            xWeight[3] = (1.0 / 6.0) * w * w * w;
            xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
            xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
            xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
            break;
        case 4L:
            /* x */
            w = x - (bspline_float_type)xIndex[2];
            w2 = w * w;
            t = (1.0 / 6.0) * w2;
            xWeight[0] = 1.0 / 2.0 - w;
            xWeight[0] *= xWeight[0];
            xWeight[0] *= (1.0 / 24.0) * xWeight[0];
            t0 = w * (t - 11.0 / 24.0);
            t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
            xWeight[1] = t1 + t0;
            xWeight[3] = t1 - t0;
            xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
            xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
            break;
        case 5L:
            /* x */
            w = x - (bspline_float_type)xIndex[2];
            w2 = w * w;
            xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
            w2 -= w;
            w4 = w2 * w2;
            w -= 1.0 / 2.0;
            t = w2 * (w2 - 3.0);
            xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
            t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
            t1 = (-1.0 / 12.0) * w * (t + 4.0);
            xWeight[2] = t0 + t1;
            xWeight[3] = t0 - t1;
            t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
            t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
            xWeight[1] = t0 + t1;
            xWeight[4] = t0 - t1;
            break;
        case 6L:
            /* x */
            w = x - (bspline_float_type)xIndex[3];
            xWeight[0] = 1.0 / 2.0 - w;
            xWeight[0] *= xWeight[0] * xWeight[0];
            xWeight[0] *= xWeight[0] / 720.0;
            xWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
                * (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w)
                    * (1.0 / 2.0 + w))))) / 120.0;
            xWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
                * (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w
                    * (-1.0 + w)))))) / 48.0;
            w2 = w * w;
            xWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
                * (21.0 / 4.0 - w2))) / 36.0;
            xWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
                * (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w
                    * (1.0 + w)))))) / 48.0;
            xWeight[6] = 1.0 / 2.0 + w;
            xWeight[6] *= xWeight[6] * xWeight[6];
            xWeight[6] *= xWeight[6] / 720.0;
            xWeight[5] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
                - xWeight[4] - xWeight[6];
            break;
        case 7L:
            /* x */
            w = x - (bspline_float_type)xIndex[3];
            xWeight[0] = 1.0 - w;
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0] * xWeight[0];
            xWeight[0] *= (1.0 - w) / 5040.0;
            w2 = w * w;
            xWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
                * (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
            xWeight[2] = (397.0 / 7.0 - w * (245.0 / 3.0 + w * (-15.0 + w
                * (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
                    * (-5.0 + w))))))) / 240.0;
            xWeight[3] = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2
                * (-4.0 + w)))) / 144.0;
            xWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
                * (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
            xWeight[5] = (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w
                * (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
            xWeight[7] = w2;
            xWeight[7] *= xWeight[7] * xWeight[7];
            xWeight[7] *= w / 5040.0;
            xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
                - xWeight[4] - xWeight[5] - xWeight[7];
            break;
        case 8L:
            /* x */
            w = x - (bspline_float_type)xIndex[4];
            xWeight[0] = 1.0 / 2.0 - w;
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0] / 40320.0;
            w2 = w * w;
            xWeight[1] = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2)))
                * (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w
                    * (-3.0 + w)))) / 5040.0;
            xWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
                * (2275.0 / 16.0 + w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w
                    * (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) / 1440.0;
            xWeight[3] = (310661.0 / 1792.0 - w * (14219.0 / 64.0 + w
                * (-199.0 / 8.0 + w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w
                    * (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) / 720.0;
            xWeight[4] = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2
                * (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
            xWeight[5] = (310661.0 / 1792.0 - w * (-14219.0 / 64.0 + w
                * (-199.0 / 8.0 + w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w
                    * (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) / 720.0;
            xWeight[7] = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2)))
                * (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w
                    * (3.0 + w)))) / 5040.0;
            xWeight[8] = 1.0 / 2.0 + w;
            xWeight[8] *= xWeight[8];
            xWeight[8] *= xWeight[8];
            xWeight[8] *= xWeight[8] / 40320.0;
            xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
                - xWeight[4] - xWeight[5] - xWeight[7] - xWeight[8];
            break;
        case 9L:
            /* x */
            w = x - (bspline_float_type)xIndex[4];
            xWeight[0] = 1.0 - w;
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0] * (1.0 - w) / 362880.0;
            xWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
                * (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 + w
                    * (24.0 + w * (-8.0 + w))))))))) / 40320.0;
            xWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 + w * (-952.0 + w
                * (938.0 / 3.0 + w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w
                    * (14.0 + w * (-7.0 + w))))))))) / 10080.0;
            xWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 + w * (66.0 + w
                * (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w
                    * (-6.0 + w))))))))) / 4320.0;
            w2 = w * w;
            xWeight[4] = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2
                * (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
            xWeight[5] = (44117.0 / 63.0 + w * (809.0 + w * (44.0 + w
                * (-868.0 / 3.0 + w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w
                    * (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
            xWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 + w * (-408.0 + w
                * (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w)
                    * (2.0 + w))))))) / 4320.0;
            xWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 + w * (118.0 + w
                * (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w
                    * (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
            xWeight[9] = w2 * w2;
            xWeight[9] *= xWeight[9] * w / 362880.0;
            xWeight[8] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
                - xWeight[4] - xWeight[5] - xWeight[6] - xWeight[7] - xWeight[9];
            break;
        default:
            GERROR_STREAM("Invalid spline degree " << SplineDegree);
        }
    }

    template <typename T, unsigned int D, typename coord_type>
    inline void hoNDBSpline<T, D, coord_type>::BSplineDiscreteFirstOrderDerivative(bspline_float_type x, unsigned int SplineDegree, bspline_float_type * weight, long long* xIndex)
    {
        unsigned int k;
        for (k = 0; k < SplineDegree; k++)
        {
            weight[k] = BSplineFirstOrderDerivative(x - xIndex[k], SplineDegree);
        }
    }

    template <typename T, unsigned int D, typename coord_type>
    inline void hoNDBSpline<T, D, coord_type>::BSplineDiscreteSecondOrderDerivative(bspline_float_type x, unsigned int SplineDegree, bspline_float_type * weight, long long* xIndex)
    {
        unsigned int k;
        for (k = 0; k < SplineDegree; k++)
        {
            weight[k] = BSplineSecondOrderDerivative(x - xIndex[k], SplineDegree);
        }
    }

    template <typename T, unsigned int D, typename coord_type>
    inline typename hoNDBSpline<T, D, coord_type>::bspline_float_type hoNDBSpline<T, D, coord_type>::BSplineFirstOrderDerivative(bspline_float_type x, unsigned int SplineDegree)
    {
        return (BSpline(x + 0.5, SplineDegree - 1) - BSpline(x - 0.5, SplineDegree - 1));
    }

    template <typename T, unsigned int D, typename coord_type>
    inline typename hoNDBSpline<T, D, coord_type>::bspline_float_type hoNDBSpline<T, D, coord_type>::BSplineSecondOrderDerivative(bspline_float_type x, unsigned int SplineDegree)
    {
        return (BSpline(x + 1, SplineDegree - 2) + BSpline(x - 1, SplineDegree - 2) - 2 * BSpline(x, SplineDegree - 2));
    }

    template <typename T, unsigned int D, typename coord_type>
    inline void hoNDBSpline<T, D, coord_type>::computeBSplineInterpolationLocationsAndWeights(size_t len, unsigned int SplineDegree, unsigned int dx, coord_type x, bspline_float_type * weight, long long* xIndex)
    {
        BSplineInterpolationLocation(x, SplineDegree, xIndex);

        if (dx == 0)
        {
            BSplineDiscrete(x, SplineDegree, weight, xIndex);
        }
        else if (dx == 1)
        {
            BSplineDiscreteFirstOrderDerivative(x, SplineDegree, weight, xIndex);
        }
        else if (dx == 2)
        {
            BSplineDiscreteSecondOrderDerivative(x, SplineDegree, weight, xIndex);
        }
        else
        {
            GERROR_STREAM("Derivative order must be 0/1/2 ... ");
            return;
        }

        BSplineInterpolationMirrorBoundaryCondition(SplineDegree, xIndex, len);
    }
}
