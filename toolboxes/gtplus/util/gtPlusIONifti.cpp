/** \file       gtPlusIONifti.cpp
    \brief      Implement the suppor for the Nifti medical image format
    \author     Hui Xue

    Ref to:
    http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
*/

#include <gtPlusIONifti.h>

// to suppor the ISMRMRD format
// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]

namespace Gadgetron { namespace gtPlus {

gtPlusIONifti::gtPlusIONifti() : BaseClass()
{
}

gtPlusIONifti::gtPlusIONifti(float px, float py) : BaseClass(px, py)
{
}

gtPlusIONifti::gtPlusIONifti(float px, float py, float pz) : BaseClass(px, py, pz)
{
}

gtPlusIONifti::gtPlusIONifti(float px, float py, float pz, float pt) : BaseClass(px, py, pz, pt)
{
}

gtPlusIONifti::gtPlusIONifti(float px, float py, float pz, float pt, float pr) : BaseClass(px, py, pz, pt, pr)
{
}

gtPlusIONifti::gtPlusIONifti(float px, float py, float pz, float pt, float pr, float ps) : BaseClass(px, py, pz, pt, pr, ps)
{
}

gtPlusIONifti::gtPlusIONifti(float px, float py, float pz, float pt, float pr, float ps, float pp) : BaseClass(px, py, pz, pt, pr, ps, pp)
{
}

gtPlusIONifti::gtPlusIONifti(float px, float py, float pz, float pt, float pr, float ps, float pp, float pq) : BaseClass(px, py, pz, pt, pr, ps, pp, pq)
{
}

void gtPlusIONifti::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus Array input/output to Nifti format -------------" << endl;
    os << "--------------------------------------------------------------------------" << endl;
}

}}
