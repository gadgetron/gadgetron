#include "mri_core.h"
#include "hoNDArray_math.h"

namespace Gadgetron{

// combine (point-wise inner product)
template <typename T> EXPORTMRICORE
void coilmap_combine(const hoNDArray<std::complex<T> > & coilmap, const hoNDArray<std::complex<T> > & data, hoNDArray<std::complex<T> > & result)
{

    //check the dimensions
    if (!(coilmap.dimensions_equal(&data))) {
        throw std::runtime_error("Coil map and data do not have the same dimensions.\n");
    }

    size_t ndims = coilmap.get_number_of_dimensions() - 1;
    if ( (ndims != 2) || (ndims != 3) ) {
        throw std::runtime_error("coilmap_combine only supports 2 or 3 dimensions.\n");
    }
    if (result.get_number_of_dimensions() != ndims) {
        throw std::runtime_error("result dimensions are not compatible with coilmap and data.\n");
    }

    size_t RO = coilmap.get_size(0);
    if (result.get_size(0) != RO) {
        throw std::runtime_error("result dimensions are not compatible with coilmap and data.\n");
    }

    size_t E1 = coilmap.get_size(1);
    if (result.get_size(1) != E1) {
        throw std::runtime_error("result dimensions are not compatible with coilmap and data.\n");
    }
    
    size_t E2 = 0;
    size_t CHA = 0;
    size_t vsize = 0;
    std::vector<size_t> vdims;
    vdims.push_back(RO);
    vdims.push_back(E1);
    if (ndims == 2) {
        CHA = coilmap.get_size(2);
        vsize = RO*E1;
    } else {
        E2 = coilmap.get_size(2);
        if (result.get_size(2) != E2) {
            throw std::runtime_error("result dimensions are not compatible with coilmap and data.\n");
        }
        CHA = coilmap.get_size(3);
        vsize = RO*E1*E2;
        vdims.push_back(E2);
    }

    //create a temporary array
    hoNDArray<std::complex<T> > temp(vdims);

    //initialize the result
    result.fill(std::complex<T>(0.0,0.0));
    
    //loop over the channels and add
    for (size_t j=0; j<CHA; j++) {
        //make to hoNDArrays that are thin-wrappers around a channel's worth of coilmap and data
        hoNDArray<std::complex<T> > cj(vdims, const_cast<std::complex<T>*>(coilmap.begin()+j*vsize));
        hoNDArray<std::complex<T> > dj(vdims, const_cast<std::complex<T>*>(data.begin()+j*vsize));
        //x*conj(y)
        multiplyConj(dj, cj, temp);
        // r = x+y
        add(result, temp, result);
    }
    
}

// scale (point-wise outer product)
template <typename T> EXPORTMRICORE
void coilmap_scale(const hoNDArray<std::complex<T> > & coilmap, const hoNDArray<std::complex<T> > & rho, hoNDArray<std::complex<T> > & result)
{

    //check the dimensions
    if (!(coilmap.dimensions_equal(&result))) {
        throw std::runtime_error("Coil map and result do not have the same dimensions.\n");
    }

    size_t ndims = coilmap.get_number_of_dimensions() - 1;
    if ( (ndims != 2) || (ndims != 3) ) {
        throw std::runtime_error("coilmap_scale only supports 2 or 3 dimensions.\n");
    }
    if (rho.get_number_of_dimensions() != ndims) {
        throw std::runtime_error("rho dimensions are not compatible with coilmap and result.\n");
    }

    size_t RO = coilmap.get_size(0);
    if (rho.get_size(0) != RO) {
        throw std::runtime_error("rho dimensions are not compatible with coilmap and result.\n");
    }

    size_t E1 = coilmap.get_size(1);
    if (rho.get_size(1) != E1) {
        throw std::runtime_error("rho dimensions are not compatible with coilmap and result.\n");
    }
    
    size_t E2 = 0;
    size_t CHA = 0;
    size_t vsize = 0;
    std::vector<size_t> vdims;
    vdims.push_back(RO);
    vdims.push_back(E1);
    if (ndims == 2) {
        CHA = coilmap.get_size(2);
        vsize = RO*E1;
    } else {
        E2 = coilmap.get_size(2);
        if (rho.get_size(2) != E2) {
            throw std::runtime_error("rho dimensions are not compatible with coilmap and result.\n");
        }
        CHA = coilmap.get_size(3);
        vsize = RO*E1*E2;
        vdims.push_back(E2);
    }

    //initialize the result
    result.fill(std::complex<T>(0.0,0.0));
    
    //loop over the channels, compute the result and stuff
    for (size_t j=0; j<CHA; j++) {
        //make to hoNDArrays that are thin-wrappers around a channel's worth of coilmap and result
        hoNDArray<std::complex<T> > cj(vdims, const_cast<std::complex<T>*>(coilmap.begin()+j*vsize));
        hoNDArray<std::complex<T> > rj(vdims, result.begin()+j*vsize);
        //rj = cj.*rho
        multiply(cj, rho, rj);            
    }

}

// norm (point-wise 2norm)
template <typename T> EXPORTMRICORE
void coilmap_norm(const hoNDArray<std::complex<T> > & coilmap, hoNDArray<T> & result)
{
    //check the dimensions (result has 1 less dimension)
    std::vector<size_t> adims;
    coilmap.get_dimensions(adims);
    std::vector<size_t> vdims(adims.begin(),adims.end()-1);

    if (!(result.dimensions_equal(&vdims))) {
        throw std::runtime_error("Coil map and result have incompatible dimensions.\n");
    }

    //initialize the result
    result.fill(0.0);

    // sum of squares
    size_t CHA = *adims.rbegin();
    size_t vsize = result.get_number_of_elements();
    for (size_t j = 0; j<CHA; j++) {
        //make an hoNDArray that is a thin-wrappers around a channel's worth of data
        hoNDArray<std::complex<T> > cj(vdims, const_cast<std::complex<T>*>(coilmap.begin()+j*vsize));
        boost::shared_ptr<hoNDArray<T> > temp = abs_square( &cj );
        add(result, *temp, result);
    }
    // square root
    sqrt_inplace(&result);
}

//instantiations
template EXPORTMRICORE void coilmap_combine(const hoNDArray<std::complex<float> > & coilmap, const hoNDArray<std::complex<float> > & data, hoNDArray<std::complex<float> > & result);
template EXPORTMRICORE void coilmap_scale(const hoNDArray<std::complex<float> > & coilmap, const hoNDArray<std::complex<float> > & rho, hoNDArray<std::complex<float> > & result);
template EXPORTMRICORE void coilmap_norm(const hoNDArray<std::complex<float> > & coilmap, hoNDArray<float> & result);

}
