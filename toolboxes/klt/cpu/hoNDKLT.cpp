
#include "hoNDKLT.h"
#include "hoMatrix.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"
#include "hoArmadillo.h"

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
    this->M_.create(N, this->output_length_, V_.begin());

    return *this;
}

template<typename T>
void hoNDKLT<T>::compute_eigen_vector(const hoNDArray<T>& data, bool remove_mean)
{

    size_t NDim = data.get_number_of_dimensions();
    size_t N = data.get_size(NDim-1);

    size_t M = data.get_number_of_elements() / N;

    hoNDArray<T> data2D;
    data2D.create(M, N, const_cast<T*>(data.begin()));

    V_.create(N, N);
    E_.create(N, 1);
    Gadgetron::clear(V_);
    Gadgetron::clear(E_);

        size_t m, n;

        // compute and subtract mean from data
        hoNDArray<T> data2DNoMean;

        arma::Mat<T> Am;
        if (remove_mean)
        {
            hoNDArray<T> dataMean(1, N);
            Gadgetron::sum_over_dimension(data2D, dataMean, 0);

            Gadgetron::scal((T)(1.0 / M), dataMean);

            data2DNoMean.create(M, N);

            for (n = 0; n < N; n++)
            {
                for (m = 0; m < M; m++)
                {
                    data2DNoMean(m, n) = data2D(m, n) - dataMean(0, n);
                }
            }

            Am = as_arma_matrix(data2DNoMean);
        }
        else
        {
            Am = as_arma_matrix(data2D);
        }

        // call svd
        arma::Mat<T> Vm = as_arma_matrix(V_);
        arma::Mat<T> Um;
        arma::Col<value_type> Sv;



        arma::svd_econ(Um, Sv, Vm, Am, 'r');


        for (n = 0; n < N; n++)
        {
            value_type v = Sv(n);
            E_(n) = v * v; // the E is eigen value, the square of singular value
        }

}

template<typename T>
void hoNDKLT<T>::prepare(const hoNDArray<T>& data, size_t dim, size_t output_length, bool remove_mean)
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
        this->compute_eigen_vector(data, remove_mean);
    }
    else if (K == 1)
    {
        std::vector<size_t> dimShrinked(dim + 1);
        memcpy(&dimShrinked[0], &dimD[0], sizeof(size_t)*(dim + 1));

        hoNDArray<T> dataS;
        dataS.create(dimShrinked, const_cast<T*>(data.begin()));

        this->compute_eigen_vector(dataS, remove_mean);
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
        Gadgetron::permute( data, dataP, dimOrder);

        this->compute_eigen_vector(dataP, remove_mean);
    }

    GADGET_CHECK_THROW(V_.get_size(0)==N);
    GADGET_CHECK_THROW(V_.get_size(1) == N);
    GADGET_CHECK_THROW(E_.get_size(0) == N);

    M_.create(N, output_length_, V_.begin());

}

template<typename T>
void hoNDKLT<T>::compute_num_kept(value_type thres)
{
    size_t N = E_.get_size(0);

    if (thres <= 0)
    {
        output_length_ = N;
    }
    else
    {
        size_t n;
        for (n = 1; n < N; n++)
        {
            if (std::abs(E_(n)) < thres*std::abs(E_(0)))
            {
                break;
            }
        }

        output_length_ = n;
    }

    GDEBUG("NUMBER OF MODES KEPT %d \n", output_length_);
}

template<typename T>
void hoNDKLT<T>::prepare(const hoNDArray<T>& data, size_t dim, value_type thres, bool remove_mean)
{
    try
    {
        this->prepare(data, dim, (size_t)0, remove_mean);
        this->compute_num_kept(thres);
        M_.create(E_.get_size(0), output_length_, V_.begin());
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::prepare(thres) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::exclude_untransformed(const hoNDArray<T>& data, size_t dim, std::vector<size_t>& untransformed, hoNDArray<T>& dataCropped)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(dim < NDim);

        size_t N = data.get_size(dim);

        size_t unN = untransformed.size();
        GADGET_CHECK_THROW(unN < N);

        size_t d;
        for (d = 0; d < unN; d++)
        {
            GADGET_CHECK_THROW(untransformed[d] < N);
        }

        // crop the data to exclude untransformed slots
        std::vector<size_t> dims;
        data.get_dimensions(dims);

        std::vector<size_t> dimCropped(dims);
        dimCropped[dim] = N - unN;

        dataCropped.create(dimCropped);

        size_t numCopySize = 1;
        for (d = 0; d < dim; d++)
        {
            numCopySize *= dims[d];
        }

        size_t num = data.get_number_of_elements() / (numCopySize*N);

        size_t n;
        for (n = 0; n < num; n++)
        {
            size_t ind = 0;
            for (d = 0; d < N; d++)
            {
                bool isUntransformed = false;
                for (size_t un = 0; un < unN; un++)
                {
                    if (untransformed[un] == d)
                    {
                        isUntransformed = true;
                        break;
                    }
                }

                if (isUntransformed)
                {
                    continue;
                }
                else
                {
                    memcpy(dataCropped.begin() + n*numCopySize*(N - unN) + ind*numCopySize, data.begin() + n*numCopySize*N + d*numCopySize, sizeof(T)*numCopySize);
                    ind++;
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::exclude_untransformed(...) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::copy_and_reset_transform(size_t N, std::vector<size_t>& untransformed)
{
    try
    {
        // adjust the eigen vector matrix
        hoNDArray<T> V;
        V.create(N, N);
        Gadgetron::clear(V);

        hoNDArray<T> E;
        E.create(N, 1);
        Gadgetron::clear(V);

        size_t unN = untransformed.size();

        size_t d;
        // set the columns for the untransformed slots
        for (d = 0; d < unN; d++)
        {
            V(untransformed[d], d) = 1;
            E(d) = E_(0);
        }

        // set the colunmns for the transformed slots
        for (d = unN; d < N; d++)
        {
            size_t ind = 0;
            for (size_t n = 0; n < N; n++)
            {
                bool isUntransformed = false;
                for (size_t un = 0; un < unN; un++)
                {
                    if (untransformed[un] == n)
                    {
                        isUntransformed = true;
                        break;
                    }
                }

                if (!isUntransformed)
                {
                    V(n, d) = V_(ind, d - unN);
                    ind++;
                }
            }

            E(d) = E_(d - unN);
        }

        V_ = V;
        E_ = E;
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::copy_and_reset_transform(...) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::prepare(const hoNDArray<T>& data, size_t dim, std::vector<size_t>& untransformed, size_t output_length, bool remove_mean)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(dim<NDim);

        size_t N = data.get_size(dim);

        size_t unN = untransformed.size();
        GADGET_CHECK_THROW(unN<N);
        if (output_length > 0)
        {
            GADGET_CHECK_THROW(output_length >= unN);
        }

        size_t d;
        for (d = 0; d < unN; d++)
        {
            GADGET_CHECK_THROW(untransformed[d] < N);
        }

        if (unN > 0)
        {
            // crop the data to exclude untransformed slots
            hoNDArray<T> dataCropped;
            this->exclude_untransformed(data, dim, untransformed, dataCropped);

            // compute the transform
            if (output_length > 0)
            {
                this->prepare(dataCropped, dim, output_length - unN, remove_mean);
            }
            else
            {
                this->prepare(dataCropped, dim, (size_t)0, remove_mean);
            }

            // adjust the eigen vector matrix
            this->copy_and_reset_transform(N, untransformed);

            output_length_ += unN;

            M_.create(N, output_length_, V_.begin());
        }
        else
        {
            this->prepare(data, dim, output_length, remove_mean);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::prepare(untransformed, thres) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::prepare(const hoNDArray<T>& data, size_t dim, std::vector<size_t>& untransformed, value_type thres, bool remove_mean)
{
    try
    {
        size_t unN = untransformed.size();

        if (unN > 0)
        {
            this->prepare(data, dim, untransformed, (size_t)(0), remove_mean);
            this->compute_num_kept(thres);
            M_.create(data.get_size(dim), output_length_, V_.begin());
        }
        else
        {
            this->prepare(data, dim, thres, remove_mean);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::prepare(untransformed, thres) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::transform(const hoNDArray<T>& in, hoNDArray<T>& out, size_t dim) const
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
            Gadgetron::permute(in, inP, dimOrder);
            hoNDArray<T> inP2D;
            inP2D.create(num, N, inP.begin());

            dimPermuted[NDim - 1] = M_.get_size(1);
            hoNDArray<T> outP;
            outP.create(dimPermuted);
            hoNDArray<T> outP2D;
            outP2D.create(num, M_.get_size(1), outP.begin());

            Gadgetron::gemm(outP2D, inP2D, false, M_, false);

            Gadgetron::permute(outP, out, dimOrder);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::transform(hoNDArray<T>& in, hoNDArray<T>& out, size_t dim) ... ");
    }
}

template<typename T>
void hoNDKLT<T>::KL_filter(const hoNDArray<T>& in, hoNDArray<T>& out, size_t dim, size_t mode_kept) const
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
        for (c = mode_kept; c<N; c++)
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
            Gadgetron::permute(in, inP, dimOrder);
            hoNDArray<T> inP2D;
            inP2D.create(num, N, inP.begin());

            hoNDArray<T> outP;
            outP.create(dimPermuted);
            hoNDArray<T> outP2D;
            outP2D.create(num, N, outP.begin());

            Gadgetron::gemm(outP2D, inP2D, false, EET, false);

            Gadgetron::permute(outP, out, dimOrder);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDKLT<T>::KL_filter(const hoNDArray<T>& in, hoNDArray<T>& out, size_t dim, size_t mode_kept) ... ");
    }
}

template<typename T>
size_t hoNDKLT<T>::transform_length() const
{
    return M_.get_size(0);
}

template<typename T>
size_t hoNDKLT<T>::output_length() const
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

        M_.create(N, output_length_, V_.begin());
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

            memcpy(M_.begin(), Mc.begin(), sizeof(T)*N*length);
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

template class hoNDKLT<float>;
template class hoNDKLT<double>;
template class hoNDKLT< std::complex<float> >;
template class hoNDKLT< std::complex<double> >;
}
