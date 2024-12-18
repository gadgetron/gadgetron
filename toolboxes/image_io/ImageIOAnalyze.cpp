/** \file       ImageIOAnalyze.cpp
    \brief      Implement the suppor for the Analzye75 medical image format
    \author     Hui Xue

    Ref to:
    http://eeg.sourceforge.net/ANALYZE75.pdf
*/

#include "ImageIOAnalyze.h"

// to support the MRD format
// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]

namespace Gadgetron {

ImageIOAnalyze::ImageIOAnalyze() : BaseClass()
{
}

ImageIOAnalyze::ImageIOAnalyze(float px, float py) : BaseClass(px, py)
{
}

ImageIOAnalyze::ImageIOAnalyze(float px, float py, float pz) : BaseClass(px, py, pz)
{
}

ImageIOAnalyze::ImageIOAnalyze(float px, float py, float pz, float pt) : BaseClass(px, py, pz, pt)
{
}

ImageIOAnalyze::ImageIOAnalyze(float px, float py, float pz, float pt, float pr) : BaseClass(px, py, pz, pt, pr)
{
}

ImageIOAnalyze::ImageIOAnalyze(float px, float py, float pz, float pt, float pr, float ps) : BaseClass(px, py, pz, pt, pr, ps)
{
}

ImageIOAnalyze::ImageIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp) : BaseClass(px, py, pz, pt, pr, ps, pp)
{
}

ImageIOAnalyze::ImageIOAnalyze(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq) : BaseClass(px, py, pz, pt, pr, ps, pp, pq)
{
}

bool ImageIOAnalyze::read_header(const std::string& filename, HeaderType& header)
{
    try
    {
        std::string filenameData = filename;
        filenameData.append(".hdr");

        ImageIOWorker ioworker(filenameData, true);

        GADGET_CHECK_RETURN_FALSE(ioworker.open());
        GADGET_CHECK_RETURN_FALSE(ioworker.read(reinterpret_cast<char*>(&header), sizeof(dsr)));
        GADGET_CHECK_RETURN_FALSE(ioworker.close());
    }
    catch(...)
    {
        GERROR_STREAM("Errors in ImageIOAnalyze::read_header(const std::string& filename, dsr& header) ... ");
        return false;
    }

    return true;
}

bool ImageIOAnalyze::write_header(const std::string& filename, const HeaderType& header)
{
    try
    {
        std::string filenameData = filename;
        filenameData.append(".hdr");

        ImageIOWorker ioworker(filenameData, false);

        GADGET_CHECK_RETURN_FALSE(ioworker.open());
        GADGET_CHECK_RETURN_FALSE(ioworker.write(reinterpret_cast<const char*>(&header), sizeof(dsr)));
        GADGET_CHECK_RETURN_FALSE(ioworker.close());
    }
    catch(...)
    {
        GERROR_STREAM("Errors in ImageIOAnalyze::write_header(const std::string& filename, const dsr& header) ... ");
        return false;
    }

    return true;
}

}
