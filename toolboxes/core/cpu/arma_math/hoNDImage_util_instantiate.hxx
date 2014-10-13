
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<float, DimImage>& img, float sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<double, DimImage>& img, double sigma[], double* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<float, DimImage>& img, double sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<double, DimImage>& img, float sigma[], double* mem);

template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<GT_Complex8, DimImage>& img, float sigma[], GT_Complex8* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<GT_Complex16, DimImage>& img, double sigma[], GT_Complex16* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<GT_Complex8, DimImage>& img, double sigma[], GT_Complex8* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<GT_Complex16, DimImage>& img, float sigma[], GT_Complex16* mem);

template EXPORTCPUCOREMATH bool gradient(const hoNDImage<float, DimImage>& x, hoNDImage<float, DimImage> gx[]);
template EXPORTCPUCOREMATH bool gradient(const hoNDImage<double, DimImage>& x, hoNDImage<double, DimImage> gx[]);
template EXPORTCPUCOREMATH bool gradient(const hoNDImage<GT_Complex8, DimImage>& x, hoNDImage<GT_Complex8, DimImage> gx[]);
template EXPORTCPUCOREMATH bool gradient(const hoNDImage<GT_Complex16, DimImage>& x, hoNDImage<GT_Complex16, DimImage> gx[]);

#ifdef USE_MKL 

template EXPORTCPUCOREMATH bool add(const hoNDImage<float, DimImage>& x, const hoNDImage<float, DimImage>& y, hoNDImage<float, DimImage>& r);
template EXPORTCPUCOREMATH bool subtract(const hoNDImage<float, DimImage>& x, const hoNDImage<float, DimImage>& y, hoNDImage<float, DimImage>& r);
template EXPORTCPUCOREMATH bool multiply(const hoNDImage<float, DimImage>& x, const hoNDImage<float, DimImage>& y, hoNDImage<float, DimImage>& r);
template EXPORTCPUCOREMATH bool divide(const hoNDImage<float, DimImage>& x, const hoNDImage<float, DimImage>& y, hoNDImage<float, DimImage>& r);
template EXPORTCPUCOREMATH bool absolute(const hoNDImage<float, DimImage>& x, hoNDImage<float, DimImage>& r);
template EXPORTCPUCOREMATH bool argument(const hoNDImage<float, DimImage>& x, hoNDImage<float, DimImage>& r);
template EXPORTCPUCOREMATH bool sqrt(const hoNDImage<float, DimImage>& x, hoNDImage<float, DimImage>& r);
template EXPORTCPUCOREMATH bool minAbsolute(const hoNDImage<float, DimImage>& x, float& r, size_t& ind);
template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDImage<float, DimImage>& x, float& r, size_t& ind);
template EXPORTCPUCOREMATH bool addEpsilon(hoNDImage<float, DimImage>& x);
template EXPORTCPUCOREMATH bool norm2(const hoNDImage<float, DimImage>& x, float& r);
template EXPORTCPUCOREMATH bool norm1(const hoNDImage<float, DimImage>& x, float& r);
template EXPORTCPUCOREMATH bool conv2(const hoNDImage<float, DimImage>& x, const hoNDImage<float, DimImage>& y, hoNDImage<float, DimImage>& z);
template EXPORTCPUCOREMATH bool conv3(const hoNDImage<float, DimImage>& x, const hoNDImage<float, DimImage>& y, hoNDImage<float, DimImage>& z);
template EXPORTCPUCOREMATH bool inv(const hoNDImage<float, DimImage>& x, hoNDImage<float, DimImage>& r);

template EXPORTCPUCOREMATH bool add(const hoNDImage<double, DimImage>& x, const hoNDImage<double, DimImage>& y, hoNDImage<double, DimImage>& r);
template EXPORTCPUCOREMATH bool subtract(const hoNDImage<double, DimImage>& x, const hoNDImage<double, DimImage>& y, hoNDImage<double, DimImage>& r);
template EXPORTCPUCOREMATH bool multiply(const hoNDImage<double, DimImage>& x, const hoNDImage<double, DimImage>& y, hoNDImage<double, DimImage>& r);
template EXPORTCPUCOREMATH bool divide(const hoNDImage<double, DimImage>& x, const hoNDImage<double, DimImage>& y, hoNDImage<double, DimImage>& r);
template EXPORTCPUCOREMATH bool absolute(const hoNDImage<double, DimImage>& x, hoNDImage<double, DimImage>& r);
template EXPORTCPUCOREMATH bool argument(const hoNDImage<double, DimImage>& x, hoNDImage<double, DimImage>& r);
template EXPORTCPUCOREMATH bool sqrt(const hoNDImage<double, DimImage>& x, hoNDImage<double, DimImage>& r);
template EXPORTCPUCOREMATH bool minAbsolute(const hoNDImage<double, DimImage>& x, double& r, size_t& ind);
template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDImage<double, DimImage>& x, double& r, size_t& ind);
template EXPORTCPUCOREMATH bool addEpsilon(hoNDImage<double, DimImage>& x);
template EXPORTCPUCOREMATH bool norm2(const hoNDImage<double, DimImage>& x, double& r);
template EXPORTCPUCOREMATH bool norm1(const hoNDImage<double, DimImage>& x, double& r);
template EXPORTCPUCOREMATH bool conv2(const hoNDImage<double, DimImage>& x, const hoNDImage<double, DimImage>& y, hoNDImage<double, DimImage>& z);
template EXPORTCPUCOREMATH bool conv3(const hoNDImage<double, DimImage>& x, const hoNDImage<double, DimImage>& y, hoNDImage<double, DimImage>& z);
template EXPORTCPUCOREMATH bool inv(const hoNDImage<double, DimImage>& x, hoNDImage<double, DimImage>& r);

template EXPORTCPUCOREMATH bool add(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, hoNDImage<GT_Complex8, DimImage>& r);
template EXPORTCPUCOREMATH bool subtract(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, hoNDImage<GT_Complex8, DimImage>& r);
template EXPORTCPUCOREMATH bool multiply(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, hoNDImage<GT_Complex8, DimImage>& r);
template EXPORTCPUCOREMATH bool divide(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, hoNDImage<GT_Complex8, DimImage>& r);
template EXPORTCPUCOREMATH bool absolute(const hoNDImage<GT_Complex8, DimImage>& x, hoNDImage<float, DimImage>& r);
template EXPORTCPUCOREMATH bool absolute(const hoNDImage<GT_Complex8, DimImage>& x, hoNDImage<GT_Complex8, DimImage>& r);
template EXPORTCPUCOREMATH bool sqrt(const hoNDImage<GT_Complex8, DimImage>& x, hoNDImage<GT_Complex8, DimImage>& r);
template EXPORTCPUCOREMATH bool minAbsolute(const hoNDImage<GT_Complex8, DimImage>& x, GT_Complex8& r, size_t& ind);
template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDImage<GT_Complex8, DimImage>& x, GT_Complex8& r, size_t& ind);
template EXPORTCPUCOREMATH bool multiplyConj(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, hoNDImage<GT_Complex8, DimImage>& r); // r = x * conj(y)
template EXPORTCPUCOREMATH bool argument(const hoNDImage<GT_Complex8, DimImage>& x, hoNDImage<float, DimImage>& r); // r = angle(x)
template EXPORTCPUCOREMATH bool conjugate(const hoNDImage<GT_Complex8, DimImage>& x, hoNDImage<GT_Complex8, DimImage>& r); // r = conj(x)
template EXPORTCPUCOREMATH bool addEpsilon(hoNDImage<GT_Complex8, DimImage>& x);
template EXPORTCPUCOREMATH bool norm2(const hoNDImage<GT_Complex8, DimImage>& x, float& r);
template EXPORTCPUCOREMATH bool norm1(const hoNDImage<GT_Complex8, DimImage>& x, float& r);
template EXPORTCPUCOREMATH bool dotc(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, GT_Complex8& r); // x'*y, x and y are N*1 vector
template EXPORTCPUCOREMATH bool conv2(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, hoNDImage<GT_Complex8, DimImage>& z);
template EXPORTCPUCOREMATH bool conv3(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, hoNDImage<GT_Complex8, DimImage>& z);
template EXPORTCPUCOREMATH bool corr2(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, hoNDImage<GT_Complex8, DimImage>& z); // x: input data [RO E1 ...], y: corr kernel [kro ke1], z: output; each 2D slice is correlated
template EXPORTCPUCOREMATH bool corr3(const hoNDImage<GT_Complex8, DimImage>& x, const hoNDImage<GT_Complex8, DimImage>& y, hoNDImage<GT_Complex8, DimImage>& z); // x: input data [RO E1 E2 ...], y: corr kernel [kro ke1 ke2], z: output; each 3D volume is correlated
template EXPORTCPUCOREMATH bool inv(const hoNDImage<GT_Complex8, DimImage>& x, hoNDImage<GT_Complex8, DimImage>& r);

template EXPORTCPUCOREMATH bool add(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, hoNDImage<GT_Complex16, DimImage>& r);
template EXPORTCPUCOREMATH bool subtract(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, hoNDImage<GT_Complex16, DimImage>& r);
template EXPORTCPUCOREMATH bool multiply(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, hoNDImage<GT_Complex16, DimImage>& r);
template EXPORTCPUCOREMATH bool divide(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, hoNDImage<GT_Complex16, DimImage>& r);
template EXPORTCPUCOREMATH bool absolute(const hoNDImage<GT_Complex16, DimImage>& x, hoNDImage<double, DimImage>& r);
template EXPORTCPUCOREMATH bool absolute(const hoNDImage<GT_Complex16, DimImage>& x, hoNDImage<GT_Complex16, DimImage>& r);
template EXPORTCPUCOREMATH bool sqrt(const hoNDImage<GT_Complex16, DimImage>& x, hoNDImage<GT_Complex16, DimImage>& r);
template EXPORTCPUCOREMATH bool minAbsolute(const hoNDImage<GT_Complex16, DimImage>& x, GT_Complex16& r, size_t& ind);
template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDImage<GT_Complex16, DimImage>& x, GT_Complex16& r, size_t& ind);
template EXPORTCPUCOREMATH bool multiplyConj(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, hoNDImage<GT_Complex16, DimImage>& r);
template EXPORTCPUCOREMATH bool argument(const hoNDImage<GT_Complex16, DimImage>& x, hoNDImage<double, DimImage>& r);
template EXPORTCPUCOREMATH bool conjugate(const hoNDImage<GT_Complex16, DimImage>& x, hoNDImage<GT_Complex16, DimImage>& r);
template EXPORTCPUCOREMATH bool addEpsilon(hoNDImage<GT_Complex16, DimImage>& x);
template EXPORTCPUCOREMATH bool norm2(const hoNDImage<GT_Complex16, DimImage>& x, double& r);
template EXPORTCPUCOREMATH bool norm1(const hoNDImage<GT_Complex16, DimImage>& x, double& r);
template EXPORTCPUCOREMATH bool dotc(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, GT_Complex16& r);
template EXPORTCPUCOREMATH bool conv2(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, hoNDImage<GT_Complex16, DimImage>& z);
template EXPORTCPUCOREMATH bool conv3(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, hoNDImage<GT_Complex16, DimImage>& z);
template EXPORTCPUCOREMATH bool corr2(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, hoNDImage<GT_Complex16, DimImage>& z);
template EXPORTCPUCOREMATH bool corr3(const hoNDImage<GT_Complex16, DimImage>& x, const hoNDImage<GT_Complex16, DimImage>& y, hoNDImage<GT_Complex16, DimImage>& z);
template EXPORTCPUCOREMATH bool inv(const hoNDImage<GT_Complex16, DimImage>& x, hoNDImage<GT_Complex16, DimImage>& r);

#endif // USE_MKL
