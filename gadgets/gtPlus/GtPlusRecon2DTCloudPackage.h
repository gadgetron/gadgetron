/** \file   GtPlusRecon2DTCloudPackage.h
    \brief  To support the dual layer GtPlus cloud, this cloud job type is defined here

            Ref to: 

            Hui Xue, Souheil Inati, Thomas Sangild Sorensen, Peter Kellman, Michael S. Hansen. 
            Distributed MRI Reconstruction using Gadgetron based Cloud Computing. 
            Magenetic Resonance in Medicine, doi: 10.1002/mrm.25213.

    \author Hui Xue
*/

#pragma once

namespace Gadgetron
{

struct EXPORTGTPLUSGADGET GtPlusRecon2DTPara
{
    size_t reconSizeRO_;
    size_t reconSizeE1_;
    size_t reconSizeE2_;

    float encodingFOV_RO_;
    float encodingFOV_E1_;
    float encodingFOV_E2_;

    float reconFOV_RO_;
    float reconFOV_E1_;
    float reconFOV_E2_;

    Gadgetron::gtPlus::ISMRMRDDIM dim_4th_;
    Gadgetron::gtPlus::ISMRMRDDIM dim_5th_;
    Gadgetron::gtPlus::ISMRMRDDIM workOrder_ShareDim_;

    bool no_acceleration_averageall_ref_;
    int no_acceleration_ref_numOfModes_;
    bool no_acceleration_same_combinationcoeff_allS_;
    int no_acceleration_whichS_combinationcoeff_;

    bool interleaved_same_combinationcoeff_allS_;
    int interleaved_whichS_combinationcoeff_;
    int interleaved_ref_numOfModes_;

    bool embedded_averageall_ref_;
    int embedded_ref_numOfModes_;
    bool embedded_fullres_coilmap_;
    bool embedded_fullres_coilmap_useHighestSignal_;
    bool embedded_same_combinationcoeff_allS_;
    int embedded_whichS_combinationcoeff_;
    bool embedded_ref_fillback_;

    bool separate_averageall_ref_;
    int separate_ref_numOfModes_;
    bool separate_fullres_coilmap_;
    bool separate_same_combinationcoeff_allS_;
    int separate_whichS_combinationcoeff_;

    bool same_coil_compression_coeff_allS_;

    bool recon_kspace_needed_;

    Gadgetron::gtPlus::gtPlusReconWorkOrderPara workOrderPara_;
};

template <typename T> 
struct GtPlusRecon2DTCloudPackage
{
    typedef typename realType<T>::Type real_value_type;

    GtPlusRecon2DTPara para;

    hoNDArray<T> kspace;
    hoNDArray<real_value_type> timeStamp;
    hoNDArray<real_value_type> physioTimeStamp;

    hoNDArray<T> ref;

    hoNDArray<T> complexIm;
    hoNDArray<T> res;

    // extra recon images
    hoNDArray<T> complexImSecond;
    // optional time stamps for the secod recon images, in the unit of seconds
    hoNDArray<real_value_type> resTimeStampSecond;
    hoNDArray<real_value_type> resPhysioTimeStampSecond;

    GtPlusRecon2DTCloudPackage();
    GtPlusRecon2DTCloudPackage(const GtPlusRecon2DTCloudPackage& pack);

    ~GtPlusRecon2DTCloudPackage();

    GtPlusRecon2DTCloudPackage<T>& operator=(const GtPlusRecon2DTCloudPackage<T>& pack);

    virtual bool serialize(char*& buf, size_t& len) const ;
    virtual bool deserialize(char* buf, size_t& len);
};

template <typename T> 
GtPlusRecon2DTCloudPackage<T>::GtPlusRecon2DTCloudPackage()
{
    kspace.clear();
    timeStamp.clear();
    physioTimeStamp.clear();
    ref.clear();
    complexIm.clear();
    complexImSecond.clear();
    resTimeStampSecond.clear();
    resPhysioTimeStampSecond.clear();
    res.clear();
}

template <typename T> 
GtPlusRecon2DTCloudPackage<T>::~GtPlusRecon2DTCloudPackage()
{

}

template <typename T> 
GtPlusRecon2DTCloudPackage<T>::GtPlusRecon2DTCloudPackage(const GtPlusRecon2DTCloudPackage& pack)
{
    para = pack.para;
    kspace = pack.kspace;
    timeStamp = pack.timeStamp;
    physioTimeStamp = pack.physioTimeStamp;
    ref = pack.ref;
    complexIm = pack.complexIm;
    complexImSecond = pack.complexImSecond;
    resTimeStampSecond = pack.resTimeStampSecond;
    resPhysioTimeStampSecond = pack.resPhysioTimeStampSecond;
    res = pack.res;
}

template <typename T> 
GtPlusRecon2DTCloudPackage<T>& GtPlusRecon2DTCloudPackage<T>::operator=(const GtPlusRecon2DTCloudPackage& pack)
{
    if ( this == &pack ) return *this;

    para = pack.para;
    kspace = pack.kspace;
    timeStamp = pack.timeStamp;
    physioTimeStamp = pack.physioTimeStamp;
    ref = pack.ref;
    complexIm = pack.complexIm;
    complexImSecond = pack.complexImSecond;
    resTimeStampSecond = pack.resTimeStampSecond;
    resPhysioTimeStampSecond = pack.resPhysioTimeStampSecond;
    res = pack.res;

    return *this;
}

template <typename T> 
bool GtPlusRecon2DTCloudPackage<T>::serialize(char*& buf, size_t& len) const 
{
    char *bufKSpace(NULL), *bufTimeStamp(NULL), *bufPhysioTimeStamp(NULL), *bufRef(NULL), *bufComplexIm(NULL), *bufRes(NULL), *bufComplexImSecond(NULL), *bufResTimeStampSecond(NULL), *bufResPhysioTimeStampSecond(NULL);
    try
    {
        if ( buf != NULL ) delete[] buf;

        // find the total len
        size_t lenKSpace, lenTimeStamp, lenPhysioTimeStamp, lenRef, lenComplexIm, lenRes, lenComplexImSecond, lenResTimeStampSecond, lenResPhyisoTimeStampSecond;

        GADGET_CHECK_THROW(kspace.serialize(bufKSpace, lenKSpace));
        GADGET_CHECK_THROW(timeStamp.serialize(bufTimeStamp, lenTimeStamp));
        GADGET_CHECK_THROW(physioTimeStamp.serialize(bufPhysioTimeStamp, lenPhysioTimeStamp));
        GADGET_CHECK_THROW(ref.serialize(bufRef, lenRef));
        GADGET_CHECK_THROW(complexIm.serialize(bufComplexIm, lenComplexIm));
        GADGET_CHECK_THROW(res.serialize(bufRes, lenRes));
        GADGET_CHECK_THROW(complexImSecond.serialize(bufComplexImSecond, lenComplexImSecond));
        GADGET_CHECK_THROW(resTimeStampSecond.serialize(bufResTimeStampSecond, lenResTimeStampSecond));
        GADGET_CHECK_THROW(resPhysioTimeStampSecond.serialize(bufResPhysioTimeStampSecond, lenResPhyisoTimeStampSecond));

        // total length
        len = sizeof(GtPlusRecon2DTPara) + lenTimeStamp + lenPhysioTimeStamp + lenKSpace + lenRef + lenComplexIm + lenRes + lenComplexImSecond + lenResTimeStampSecond + lenResPhyisoTimeStampSecond;

        buf = new char[len];
        GADGET_CHECK_RETURN_FALSE( buf != NULL );

        size_t offset = 0, currLen=0;

        currLen = sizeof(GtPlusRecon2DTPara);
        memcpy(buf+offset, &para, currLen);
        offset += currLen;

        currLen = lenKSpace;
        memcpy(buf+offset, bufKSpace, currLen);
        offset += currLen;
        delete [] bufKSpace;

        currLen = lenTimeStamp;
        memcpy(buf+offset, bufTimeStamp, currLen);
        offset += currLen;
        delete [] bufTimeStamp;

        currLen = lenPhysioTimeStamp;
        memcpy(buf+offset, bufPhysioTimeStamp, currLen);
        offset += currLen;
        delete [] bufPhysioTimeStamp;

        currLen = lenRef;
        memcpy(buf+offset, bufRef, currLen);
        offset += currLen;
        delete [] bufRef;

        currLen = lenComplexIm;
        memcpy(buf+offset, bufComplexIm, currLen);
        offset += currLen;
        delete [] bufComplexIm;

        currLen = lenRes;
        memcpy(buf+offset, bufRes, currLen);
        offset += currLen;
        delete [] bufRes;

        currLen = lenComplexImSecond;
        memcpy(buf+offset, bufComplexImSecond, currLen);
        offset += currLen;
        delete [] bufComplexImSecond;

        currLen = lenResTimeStampSecond;
        memcpy(buf+offset, bufResTimeStampSecond, currLen);
        offset += currLen;
        delete [] bufResTimeStampSecond;

        currLen = lenResPhyisoTimeStampSecond;
        memcpy(buf+offset, bufResPhysioTimeStampSecond, currLen);
        offset += currLen;
        delete [] bufResPhysioTimeStampSecond;
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors happened in GtPlusRecon2DTCloudPackage<T>::serialize(...) ... ");

        if ( bufKSpace != NULL ) delete [] bufKSpace;
        if ( bufTimeStamp != NULL ) delete [] bufTimeStamp;
        if ( bufPhysioTimeStamp != NULL ) delete [] bufPhysioTimeStamp;
        if ( bufRef != NULL ) delete [] bufRef;
        if ( bufComplexIm != NULL ) delete [] bufComplexIm;
        if ( bufRes != NULL ) delete [] bufRes;
        if ( bufComplexImSecond != NULL ) delete [] bufComplexImSecond;
        if ( bufResTimeStampSecond != NULL ) delete [] bufResTimeStampSecond;
        if ( bufResPhysioTimeStampSecond != NULL ) delete [] bufResPhysioTimeStampSecond;

        return false;
    }

    return true;
}

template <typename T> 
bool GtPlusRecon2DTCloudPackage<T>::deserialize(char* buf, size_t& len)
{
    try
    {
        memcpy(&para, buf, sizeof(GtPlusRecon2DTPara));

        size_t offset(sizeof(GtPlusRecon2DTPara)), currLen=0;

        GADGET_CHECK_RETURN_FALSE(kspace.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(timeStamp.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(physioTimeStamp.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(ref.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(complexIm.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(res.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(complexImSecond.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(resTimeStampSecond.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(resPhysioTimeStampSecond.deserialize(buf+offset, currLen));
        offset += currLen;

        // total length
        len = offset;
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors happended in GtPlusRecon2DTCloudPackage<T>::deserialize(...) ...");
        return false;
    }

    return true;
}

typedef GtPlusRecon2DTCloudPackage< std::complex<float> > GtPlusRecon2DTCloudPackageCPFL;

}
