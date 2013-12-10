/** \file       gtPlusIOAnalyze.cpp
    \brief      Implement the suppor for the Analzye75 medical image format
    \author     Hui Xue

    Ref to:
    http://eeg.sourceforge.net/ANALYZE75.pdf
*/

#include <gtPlusIOAnalyze.h>

// to suppor the ISMRMRD format
// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]

namespace Gadgetron { namespace gtPlus {

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py)
{
    pixelSize_.resize(2);
    pixelSize_[0] = px;
    pixelSize_[1] = py;
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz)
{
    pixelSize_.resize(3);
    pixelSize_[0] = px;
    pixelSize_[1] = py;
    pixelSize_[2] = pz;
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt)
{
    pixelSize_.resize(4);
    pixelSize_[0] = px;
    pixelSize_[1] = py;
    pixelSize_[2] = pz;
    pixelSize_[3] = pt;
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr)
{
    pixelSize_.resize(5);
    pixelSize_[0] = px;
    pixelSize_[1] = py;
    pixelSize_[2] = pz;
    pixelSize_[3] = pt;
    pixelSize_[4] = pr;
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps)
{
    pixelSize_.resize(6);
    pixelSize_[0] = px;
    pixelSize_[1] = py;
    pixelSize_[2] = pz;
    pixelSize_[3] = pt;
    pixelSize_[4] = pr;
    pixelSize_[5] = ps;
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp)
{
    pixelSize_.resize(7);
    pixelSize_[0] = px;
    pixelSize_[1] = py;
    pixelSize_[2] = pz;
    pixelSize_[3] = pt;
    pixelSize_[4] = pr;
    pixelSize_[5] = ps;
    pixelSize_[6] = pp;
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq)
{
    pixelSize_.resize(8);
    pixelSize_[0] = px;
    pixelSize_[1] = py;
    pixelSize_[2] = pz;
    pixelSize_[3] = pt;
    pixelSize_[4] = pr;
    pixelSize_[5] = ps;
    pixelSize_[6] = pp;
    pixelSize_[7] = pq;
}

void gtPlusIOAnalyze::setPixelSize(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq)
{
    pixelSize_.resize(8);
    pixelSize_[0] = px;
    pixelSize_[1] = py;
    pixelSize_[2] = pz;
    pixelSize_[3] = pt;
    pixelSize_[4] = pr;
    pixelSize_[5] = ps;
    pixelSize_[6] = pp;
    pixelSize_[7] = pq;
}

void gtPlusIOAnalyze::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus Array input/output to Analyze75 format -------------" << endl;
    os << "--------------------------------------------------------------------------" << endl;
}

std::string gtPlusIOAnalyze::getRTTIFromAnalyzeDataType(AnalyzeDataType aDT)
{
    std::string rttiID;

    switch (aDT)
    {
    case DT_UNSIGNED_CHAR :
        rttiID = typeid(unsigned char).name();
        break;

    case DT_SIGNED_SHORT :
        rttiID = typeid(short).name();
        break;

    case DT_UNSIGNED_SHORT :
        rttiID = typeid(unsigned short).name();
        break;

    case DT_SIGNED_INT :
        rttiID = typeid(int).name();
        break;

    case DT_UNSIGNED_INT :
        rttiID = typeid(size_t).name();
        break;

    case DT_FLOAT :
        rttiID = typeid(float).name();
        break;

    case DT_DOUBLE :
        rttiID = typeid(double).name();
        break;

    case DT_COMPLEX :
        rttiID = typeid(GT_Complex8).name();
        break;

    case DT_DOUBLECOMPLEX :
        rttiID = typeid(GT_Complex16).name();
        break;

    default:
        rttiID = "UNKOWN TYPE";
    }

    return rttiID;
}

AnalyzeDataType gtPlusIOAnalyze::getAnalyzeDataTypeFromRTTI(const std::string& name)
{
    AnalyzeDataType analyzeDT = DT_ANA_UNKNOWN;

    if ( name == typeid(unsigned char).name() )
    {
        analyzeDT = DT_UNSIGNED_CHAR;
    }

    if ( name == typeid(short).name() )
    {
        analyzeDT = DT_SIGNED_SHORT;
    }

    if ( name == typeid(unsigned short).name() )
    {
        analyzeDT = DT_UNSIGNED_SHORT;
    }

    if ( name == typeid(int).name() )
    {
        analyzeDT = DT_SIGNED_INT;
    }

    if ( name == typeid(size_t).name() )
    {
        analyzeDT = DT_UNSIGNED_INT;
    }

    if ( name == typeid(float).name() )
    {
        analyzeDT = DT_FLOAT;
    }

    if ( name == typeid(double).name() )
    {
        analyzeDT = DT_DOUBLE;
    }

    if ( name == typeid(GT_Complex8).name() )
    {
        analyzeDT = DT_COMPLEX;
    }

    if ( name == typeid(GT_Complex16).name() )
    {
        analyzeDT = DT_DOUBLECOMPLEX;
    }

    return analyzeDT;
}

bool gtPlusIOAnalyze::readAnalyzeHeader(const std::string& filename, dsr& header)
{
    try
    {
        std::string filenameData = filename;
        filenameData.append(".hdr");

        gtPlusIOWorker ioworker(filenameData, true);

        GADGET_CHECK_RETURN_FALSE(ioworker.open());
        GADGET_CHECK_RETURN_FALSE(ioworker.read(reinterpret_cast<char*>(&header), sizeof(dsr)));
        GADGET_CHECK_RETURN_FALSE(ioworker.close());
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::readAnalyzeHeader(const std::string& filename, dsr& header) ... ");
        return false;
    }

    return true;
}

bool gtPlusIOAnalyze::writeAnalyzeHeader(const std::string& filename, const dsr& header)
{
    try
    {
        std::string filenameData = filename;
        filenameData.append(".hdr");

        gtPlusIOWorker ioworker(filenameData, false);

        GADGET_CHECK_RETURN_FALSE(ioworker.open());
        GADGET_CHECK_RETURN_FALSE(ioworker.write(reinterpret_cast<const char*>(&header), sizeof(dsr)));
        GADGET_CHECK_RETURN_FALSE(ioworker.close());
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOAnalyze::writeAnalyzeHeader(const std::string& filename, const dsr& header) ... ");
        return false;
    }

    return true;
}

}}
