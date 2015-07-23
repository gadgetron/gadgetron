
#include "hoNDKLT.h"
#include "hoMatrix.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"

namespace Gadgetron{

template<typename T> 
hoNDKLT<T>::hoNDKLT()
{
}

template<typename T>
hoNDKLT<T>::hoNDKLT(const hoNDArray<T>& data, size_t dim, size_t output_length)
{
    this->prepare(data, dim, output_length);
}

template<typename T>
hoNDKLT<T>::hoNDKLT(const hoNDArray<T>& data, size_t dim, value_type thres)
{
    this->prepare(data, dim, thres);
}

template<typename T>
hoNDKLT<T>::hoNDKLT(const Self& v)
{
    *this = v;
}

template<typename T>
hoNDKLT<T>::~hoNDKLT()
{
}

template<typename T>
hoNDKLT<T>& hoNDKLT<T>::operator=(const Self& v)
{
    if (&v == this) return *this;

    this->V_ = v.V_;
    this->E_ = v.E_;
    this->output_length_ = v.output_length_;

    size_t N = this->V_.get_size(0);
    this->M_.create(N, this->output_length_, V_.begin() + (N - output_length_)*N);

    return *this;
}

template<typename T>
void hoNDKLT<T>::compute_eigen_vector(const hoNDArray<T>& data, hoNDArray<T>& V, hoNDArray<T>& E, bool remove_mean)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        size_t N = data.get_size(NDim-1);

        size_t M = data.get_number_of_elements() / N;

        hoNDArray<T> data2D;
        data2D.create(M, N, const_cast<T*>(data.begin()));

        V.create(N, N);
        E.create(N, 1);
        Gadgetron::clear(V);
        Gadgetron::clear(E);

        // compute data'*data
        char uplo = 'L';
        bool isAHA = true;
        Gadgetron::herk(V, data2D, uplo, isAHA);

        hoMatrix<T> m;
        m.createMatrix(N, N, V.begin());
        m.copyLowerTriToUpper();

        if (remove_mean)
        {
            // compute and subtract mean for every mode
            hoNDArray<T> mean(N, 1);
            hoNDArray<T> dataMean(1, N, mean.begin());
            Gadgetron::sum_over_dimension(data2D, dataMean, 0);

            Gadgetron::scal((T)(1.0 / M), mean);

            hoNDArray<T> MMH(N, N);
            Gadgetron::clear(MMH);

            hoNDArray<T> meanCopy(mean);
            Gadgetron::gemm(MMH, meanCopy, false, mean, true);
            Gadgetron::subtract(V, MMH, V);
        }

        // compute eigen vector and values
        Gadgetron::heev(V, E);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::compute_eigen_vector(output_length) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::prepare(const hoNDArray<T>& data, size_t dim, size_t output_length, bool remove_mean)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(dim<NDim);

        size_t N = data.get_size(dim);

        if (output_length > 0 && output_length <= N)
        {
            output_length_ = output_length;
        }
        else
        {
            output_length_ = N;
        }

        std::vector<size_t> dimD;
        data.get_dimensions(dimD);

        size_t K = 1;
        for (size_t n = dim + 1; n < NDim; n++) K *= dimD[n];

        if (dim == NDim - 1)
        {
            this->compute_eigen_vector(data, V_, E_, remove_mean);
        }
        else if (K == 1)
        {
            std::vector<size_t> dimShrinked(dim + 1);
            memcpy(&dimShrinked[0], &dimD[0], sizeof(size_t)*(dim + 1));

            hoNDArray<T> dataS;
            dataS.create(dimShrinked, const_cast<T*>(data.begin()));

            this->compute_eigen_vector(dataS, V_, E_, remove_mean);
        }
        else
        {
            std::vector<size_t> dimOrder(NDim), dimPermuted(dimD);

            size_t l;
            for (l = 0; l<NDim; l++)
            {
                dimOrder[l] = l;
                dimPermuted[l] = dimD[l];
            }

            dimOrder[dim] = NDim - 1;
            dimOrder[NDim - 1] = dim;

            dimPermuted[dim] = dimD[NDim - 1];
            dimPermuted[NDim - 1] = dimD[dim];

            hoNDArray<T> dataP;
            dataP.create(dimPermuted);
            Gadgetron::permute( const_cast<hoNDArray<T>* >(&data), &dataP, &dimOrder);

            this->compute_eigen_vector(dataP, V_, E_, remove_mean);
        }

        GADGET_CHECK_THROW(V_.get_size(0)==N);
        GADGET_CHECK_THROW(V_.get_size(1) == N);
        GADGET_CHECK_THROW(E_.get_size(0) == N);

        M_.create(N, output_length_, V_.begin() + (N - output_length_)*N);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::prepare(output_length) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::prepare(const hoNDArray<T>& data, size_t dim, value_type thres, bool remove_mean)
{
    try
    {
        this->prepare(data, dim, (size_t)0, remove_mean);

        size_t N = E_.get_size(0);

        if (thres <= 0)
        {
            output_length_ = N;
        }
        else
        {
            long long n;
            for (n = N - 2; n >= 0; n--)
            {
                if (std::abs(E_(n)) < thres*std::abs(E_(N - 1)))
                {
                    break;
                }
            }

            output_length_ = N - n - 1;
        }

        M_.create(N, output_length_, V_.begin() + (N - output_length_)*N);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::prepare(thres) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::prepare(const hoNDArray<T>& V, const hoNDArray<T>& E, size_t output_length)
{
    try
    {
        GADGET_CHECK_THROW(V.get_size(0) == V.get_size(1));
        GADGET_CHECK_THROW(V.get_size(0) == E.get_size(0));

        V_ = V;
        E_ = E;

        size_t N = V_.get_size(0);

        if (output_length > 0 && output_length <= N)
        {
            output_length_ = output_length;
        }
        else
        {
            output_length_ = N;
        }

        M_.create(N, output_length_, V_.begin() + (N - output_length_)*N);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::prepare(V) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::prepare(const hoNDArray<T>& M)
{
    try
    {
        M_ = M;
        output_length_ = M_.get_size(1);

        size_t N = M_.get_size(0);
        GADGET_CHECK_THROW(N >= output_length_);

        V_.create(N, N);
        Gadgetron::clear(V_);

        memcpy(V_.begin() + (N - output_length_)*N, M_.begin(), M_.get_number_of_bytes());
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::prepare(M) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::transform(const hoNDArray<T>& in, hoNDArray<T>& out, size_t dim)
{
    try
    {
        size_t NDim = in.get_number_of_dimensions();
        GADGET_CHECK_THROW(dim<NDim);

        GADGET_CHECK_THROW(in.get_size(dim)==M_.get_size(0));

        std::vector<size_t> dim_in;
        in.get_dimensions(dim_in);

        std::vector<size_t> dim_out;
        dim_out = dim_in;
        dim_out[dim] = M_.get_size(1);

        out.create(dim_out);

        size_t N = M_.get_size(0);
        size_t num = in.get_number_of_elements() / N;

        size_t K = 1;
        for (size_t n = dim + 1; n < NDim; n++) K *= dim_in[n];

        if ((dim == NDim - 1) || (K == 1))
        {
            hoNDArray<T> in2D;
            in2D.create(num, N, const_cast<T*>(in.begin()));

            hoNDArray<T> out2D;
            out2D.create(num, M_.get_size(1), out.begin());

            Gadgetron::gemm(out2D, in2D, false, M_, false);
        }
        else
        {
            std::vector<size_t> dimOrder(NDim), dimPermuted(dim_in);

            size_t l;
            for (l = 0; l<NDim; l++)
            {
                dimOrder[l] = l;
                dimPermuted[l] = dim_in[l];
            }

            dimOrder[dim] = NDim-1;
            dimOrder[NDim - 1] = dim;

            dimPermuted[dim] = dim_in[NDim - 1];
            dimPermuted[NDim - 1] = dim_in[dim];

            hoNDArray<T> inP;
            inP.create(dimPermuted);
            Gadgetron::permute(const_cast< hoNDArray<T>* >(&in), &inP, &dimOrder);
            hoNDArray<T> inP2D;
            inP2D.create(num, N, inP.begin());

            dimPermuted[NDim - 1] = M_.get_size(1);
            hoNDArray<T> outP;
            outP.create(dimPermuted);
            hoNDArray<T> outP2D;
            outP2D.create(num, M_.get_size(1), outP.begin());

            Gadgetron::gemm(outP2D, inP2D, false, M_, false);

            Gadgetron::permute(&outP, &out, &dimOrder);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::transform(hoNDArray<T>& in, hoNDArray<T>& out, size_t dim) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::KL_filter(const hoNDArray<T>& in, hoNDArray<T>& out, size_t dim, size_t mode_kept)
{
    try
    {
        size_t NDim = in.get_number_of_dimensions();
        GADGET_CHECK_THROW(dim<NDim);

        GADGET_CHECK_THROW(in.get_size(dim) == V_.get_size(0));

        std::vector<size_t> dim_in;
        in.get_dimensions(dim_in);

        out.create(dim_in);

        size_t N = V_.get_size(0);
        size_t num = in.get_number_of_elements() / N;

        hoMatrix<T> E(N, N);
        memcpy(E.begin(), V_.begin(), sizeof(T)*N*N);

        size_t r, c;
        for (c = 0; c<N - mode_kept + 1; c++)
        {
            for (r = 0; r<N; r++)
            {
                E(r, c) = T(0);
            }
        }

        hoMatrix<T> ET;
        ET.createMatrix(N, N);
        memcpy(ET.begin(), V_.begin(), sizeof(T)*N*N);

        Gadgetron::conjugatetrans(E, ET);

        hoMatrix<T> EET(N, N);
        Gadgetron::clear(EET);
        Gadgetron::gemm(EET, E, false, ET, false);

        size_t K = 1;
        for (size_t n = dim + 1; n < NDim; n++) K *= dim_in[n];

        if ((dim == NDim - 1) || (K==1))
        {
            hoNDArray<T> in2D;
            in2D.create(num, N, const_cast<T*>(in.begin()));

            hoNDArray<T> out2D;
            out2D.create(num, N, out.begin());

            Gadgetron::gemm(out2D, in2D, false, EET, false);
        }
        else
        {
            std::vector<size_t> dimOrder(NDim), dimPermuted(dim_in);

            size_t l;
            for (l = 0; l<NDim; l++)
            {
                dimOrder[l] = l;
                dimPermuted[l] = dim_in[l];
            }

            dimOrder[dim] = NDim - 1;
            dimOrder[NDim - 1] = dim;

            dimPermuted[dim] = dim_in[NDim - 1];
            dimPermuted[NDim - 1] = dim_in[dim];

            hoNDArray<T> inP;
            inP.create(dimPermuted);
            Gadgetron::permute(const_cast< hoNDArray<T>* >(&in), &inP, &dimOrder);
            hoNDArray<T> inP2D;
            inP2D.create(num, N, inP.begin());

            hoNDArray<T> outP;
            outP.create(dimPermuted);
            hoNDArray<T> outP2D;
            outP2D.create(num, N, outP.begin());

            Gadgetron::gemm(outP2D, inP2D, false, EET, false);

            Gadgetron::permute(&outP, &out, &dimOrder);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::KL_filter(const hoNDArray<T>& in, hoNDArray<T>& out, size_t dim, size_t mode_kept) ... ");
    }
}

template<typename T>
size_t hoNDKLT<T>::transform_length()
{
    return M_.get_size(0);
}

template<typename T>
size_t hoNDKLT<T>::output_length()
{
    GADGET_CHECK_THROW(M_.get_size(1) == output_length_);
    return output_length_;
}

template<typename T>
void hoNDKLT<T>::output_length(size_t length)
{
    if (M_.get_size(0) == V_.get_size(0))
    {
        GADGET_CHECK_THROW(V_.get_size(0) == V_.get_size(1));

        size_t N = M_.get_size(0);

        if (length > 0 && length <= N)
        {
            output_length_ = length;
        }
        else
        {
            output_length_ = N;
        }

        M_.create(N, output_length_, V_.begin() + (N - output_length_)*N);
    }
    else
    {
        size_t N = M_.get_size(0);
        size_t M = M_.get_size(1);

        output_length_ = (length <= M) ? length : M;

        if (output_length_<M)
        {
            hoNDArray<T> Mc(M_);
            M_.create(N, length);

            memcpy(M_.begin(), Mc.begin() + N*(M-output_length_), sizeof(T)*N*length);
        }
    }
}

template<typename T>
void hoNDKLT<T>::KL_transformation(hoNDArray<T>& M) const
{
    M = M_;
}

template<typename T>
void hoNDKLT<T>::eigen_vector(hoNDArray<T>& V) const
{
    V = V_;
}

template<typename T>
void hoNDKLT<T>::eigen_value(hoNDArray<T>& E) const
{
    E = E_;
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCPUKLT hoNDKLT<float>;
template class EXPORTCPUKLT hoNDKLT<double>;
template class EXPORTCPUKLT hoNDKLT< std::complex<float> >;
template class EXPORTCPUKLT hoNDKLT< std::complex<double> >;
}
