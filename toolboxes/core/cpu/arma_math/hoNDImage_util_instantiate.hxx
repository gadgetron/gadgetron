
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
