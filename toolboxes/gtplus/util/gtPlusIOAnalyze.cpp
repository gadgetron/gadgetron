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

gtPlusIOAnalyze::gtPlusIOAnalyze() : BaseClass()
{
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py) : BaseClass(px, py)
{
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz) : BaseClass(px, py, pz)
{
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt) : BaseClass(px, py, pz, pt)
{
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr) : BaseClass(px, py, pz, pt, pr)
{
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps) : BaseClass(px, py, pz, pt, pr, ps)
{
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp) : BaseClass(px, py, pz, pt, pr, ps, pp)
{
}

gtPlusIOAnalyze::gtPlusIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq) : BaseClass(px, py, pz, pt, pr, ps, pp, pq)
{
}

void gtPlusIOAnalyze::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus Array input/output to Analyze75 format -------------" << endl;
    os << "--------------------------------------------------------------------------" << endl;
}

bool gtPlusIOAnalyze::readHeader(const std::string& filename, HeaderType& header)
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
        GERROR_STREAM("Errors in gtPlusIOAnalyze::readHeader(const std::string& filename, dsr& header) ... ");
        return false;
    }

    return true;
}

bool gtPlusIOAnalyze::writeHeader(const std::string& filename, const HeaderType& header)
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
        GERROR_STREAM("Errors in gtPlusIOAnalyze::writeHeader(const std::string& filename, const dsr& header) ... ");
        return false;
    }

    return true;
}

}}
