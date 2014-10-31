
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<float, DimImage>& img, float sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<double, DimImage>& img, double sigma[], double* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<float, DimImage>& img, double sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage<double, DimImage>& img, float sigma[], double* mem);

template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage< std::complex<float> , DimImage>& img, float sigma[],  std::complex<float> * mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage< std::complex<double> , DimImage>& img, double sigma[],  std::complex<double> * mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage< std::complex<float> , DimImage>& img, double sigma[],  std::complex<float> * mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDImage< std::complex<double> , DimImage>& img, float sigma[],  std::complex<double> * mem);

template EXPORTCPUCOREMATH bool gradient(const hoNDImage<float, DimImage>& x, hoNDImage<float, DimImage> gx[]);
template EXPORTCPUCOREMATH bool gradient(const hoNDImage<double, DimImage>& x, hoNDImage<double, DimImage> gx[]);
template EXPORTCPUCOREMATH bool gradient(const hoNDImage< std::complex<float> , DimImage>& x, hoNDImage< std::complex<float> , DimImage> gx[]);
template EXPORTCPUCOREMATH bool gradient(const hoNDImage< std::complex<double> , DimImage>& x, hoNDImage< std::complex<double> , DimImage> gx[]);
