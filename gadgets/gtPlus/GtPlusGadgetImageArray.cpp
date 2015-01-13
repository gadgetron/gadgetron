
#include "GtPlusGadgetImageArray.h"

namespace Gadgetron
{

GtPlusGadgetImageExt::GtPlusGadgetImageExt() : ISMRMRD::ImageHeader()
{
    version = 0;
    flags = 0;
    measurement_uid = 0;

    matrix_size[0] = 0; matrix_size[1] = 0; matrix_size[2] = 0;
    field_of_view[0] = 0; field_of_view[1] = 0; field_of_view[2] = 0;
    channels = 0;
    memset(position, 0, sizeof(float)*ISMRMRD::ISMRMRD_POSITION_LENGTH);
    memset(read_dir, 0, sizeof(float)*ISMRMRD::ISMRMRD_POSITION_LENGTH);
    memset(phase_dir, 0, sizeof(float)*ISMRMRD::ISMRMRD_POSITION_LENGTH);
    memset(slice_dir, 0, sizeof(float)*ISMRMRD::ISMRMRD_POSITION_LENGTH);
    memset(patient_table_position, 0, sizeof(float)*ISMRMRD::ISMRMRD_POSITION_LENGTH);

    average = 0;
    slice = 0;
    contrast = 0;
    phase = 0;
    repetition = 0;
    set = 0;
    acquisition_time_stamp = 0;
    memset(physiology_time_stamp, 0, sizeof(uint32_t)*ISMRMRD::ISMRMRD_PHYS_STAMPS);

    data_type = 0;
    image_type = 0;
    image_index = 0;
    image_series_index = 0;

    memset(user_int, 0, sizeof(int32_t)*ISMRMRD::ISMRMRD_USER_INTS);
    memset(user_float, 0, sizeof(float)*ISMRMRD::ISMRMRD_USER_FLOATS);

    time_stamps.clear();
    pmu_time_stamps.clear();
}

GtPlusGadgetImageExt::~GtPlusGadgetImageExt()
{
}

void GtPlusGadgetImageExt::set_matrix_size(size_t index, size_t size)
{
    if (index < 3) 
    {
        matrix_size[index] = (uint16_t)size;
    }

    if ( index == 1 )
    {
        time_stamps.clear();
        time_stamps.resize(matrix_size[1], -1);
        pmu_time_stamps.clear();
        pmu_time_stamps.resize(matrix_size[1], -1);
    }
}

void GtPlusGadgetImageExt::copy(GtPlusGadgetImageExt& aMessageImage)
{
    version = aMessageImage.version;
    flags = aMessageImage.flags;
    measurement_uid = aMessageImage.measurement_uid;

    matrix_size[0] = aMessageImage.matrix_size[0];
    matrix_size[1] = aMessageImage.matrix_size[1];
    matrix_size[2] = aMessageImage.matrix_size[2];

    field_of_view[0] = aMessageImage.field_of_view[0];
    field_of_view[1] = aMessageImage.field_of_view[1];
    field_of_view[2] = aMessageImage.field_of_view[2];

    channels = aMessageImage.channels;

    memcpy(position, aMessageImage.position, sizeof(float)*ISMRMRD::ISMRMRD_POSITION_LENGTH);
    memcpy(read_dir, aMessageImage.read_dir, sizeof(float)*ISMRMRD::ISMRMRD_DIRECTION_LENGTH);
    memcpy(phase_dir, aMessageImage.phase_dir, sizeof(float)*ISMRMRD::ISMRMRD_DIRECTION_LENGTH);
    memcpy(slice_dir, aMessageImage.slice_dir, sizeof(float)*ISMRMRD::ISMRMRD_DIRECTION_LENGTH);
    memcpy(patient_table_position, aMessageImage.patient_table_position, sizeof(float)*ISMRMRD::ISMRMRD_POSITION_LENGTH);

    average = aMessageImage.average;
    slice = aMessageImage.slice;
    contrast = aMessageImage.contrast;
    phase = aMessageImage.phase;
    repetition = aMessageImage.repetition;
    set = aMessageImage.set;

    acquisition_time_stamp = aMessageImage.acquisition_time_stamp;

    memcpy(physiology_time_stamp, aMessageImage.physiology_time_stamp, sizeof(uint32_t)*ISMRMRD::ISMRMRD_PHYS_STAMPS);

    data_type = aMessageImage.data_type;
    image_type = aMessageImage.image_type;
    image_index = aMessageImage.image_index;
    image_series_index = aMessageImage.image_series_index;

    memcpy(user_int, aMessageImage.user_int, sizeof(int32_t)*ISMRMRD::ISMRMRD_USER_INTS);
    memcpy(user_float, aMessageImage.user_float, sizeof(float)*ISMRMRD::ISMRMRD_USER_FLOATS);

    time_stamps = aMessageImage.time_stamps;
    pmu_time_stamps = aMessageImage.pmu_time_stamps;
}

void GtPlusGadgetImageExt::recomputeHeader(const GtPlusGadgetImageExt& aMessageImage, double weight)
{
    size_t ii;
    for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
    {
        position[ii] = (float)((position[ii]*weight) + (1.0-weight)*aMessageImage.position[ii]);
        patient_table_position[ii] = (float)((patient_table_position[ii]*weight) + (1.0-weight)*aMessageImage.patient_table_position[ii]);
    }

    acquisition_time_stamp = (uint32_t)((acquisition_time_stamp*weight) + (1.0-weight)*aMessageImage.acquisition_time_stamp + 0.5);

    for ( ii=0; ii<ISMRMRD::ISMRMRD_PHYS_STAMPS; ii++ )
    {
        physiology_time_stamp[ii] = (uint32_t)((physiology_time_stamp[ii]*weight) + (1.0-weight)*aMessageImage.physiology_time_stamp[ii] + 0.5);
    }
}

void GtPlusGadgetImageExt::dump()
{
    using namespace std;

    cout << "GtPlusGadgetImageExt" << endl;
    cout << "----------------------------------------------------------" << endl;
    cout << "version            : " << version << endl;
    cout << "flags              : " << flags << endl;
    cout << "measurement_uid    : " << measurement_uid << endl;
    cout << "matrix_size[3]     : " << matrix_size[0] << " " << matrix_size[1] << " " << matrix_size[2] << endl;
    cout << "field_of_view[3]   : " << field_of_view[0] << " " << field_of_view[1] << " " << field_of_view[2] << endl;
    cout << "channels           : " << channels << endl;

    size_t ii;

    cout << "position[ISMRMRD::ISMRMRD_POSITION_LENGTH]      : ";
    for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
    {
        cout << position[ii] << " ";
    }
    cout << endl;

    cout << "read_dir[ISMRMRD::ISMRMRD_POSITION_LENGTH]      : ";
    for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
    {
        cout << read_dir[ii] << " ";
    }
    cout << endl;

    cout << "phase_dir[ISMRMRD::ISMRMRD_POSITION_LENGTH]      : ";
    for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
    {
        cout << phase_dir[ii] << " ";
    }
    cout << endl;

    cout << "slice_dir[ISMRMRD::ISMRMRD_POSITION_LENGTH]      : ";
    for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
    {
        cout << slice_dir[ii] << " ";
    }
    cout << endl;

    cout << "patient_table_position[ISMRMRD::ISMRMRD_POSITION_LENGTH]      : ";
    for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
    {
        cout << patient_table_position[ii] << " ";
    }
    cout << endl;

    cout << "average            : " << average << endl;
    cout << "slice              : " << slice << endl;
    cout << "contrast           : " << contrast << endl;
    cout << "phase              : " << phase << endl;
    cout << "repetition         : " << repetition << endl;
    cout << "set                : " << set << endl;
    cout << "acquisition_time_stamp : " << acquisition_time_stamp << endl;

    cout << "physiology_time_stamp[ISMRMRD::ISMRMRD_PHYS_STAMPS] : ";
    for ( ii=0; ii<ISMRMRD::ISMRMRD_PHYS_STAMPS; ii++ )
    {
        cout << physiology_time_stamp[ii] << " ";
    }
    cout << endl;

    cout << "data_type          : " << data_type << endl;
    cout << "image_type         : " << image_type << endl;
    cout << "image_index        : " << image_index << endl;
    cout << "image_series_index : " << image_series_index << endl;

    cout << "user_int[ISMRMRD::ISMRMRD_USER_INTS]        : ";
    for ( ii=0; ii<ISMRMRD::ISMRMRD_USER_INTS; ii++ )
    {
        cout << user_int[ii] << " ";
    }
    cout << endl;

    cout << "user_float[ISMRMRD::ISMRMRD_USER_FLOATS]    : ";
    for ( ii=0; ii<ISMRMRD::ISMRMRD_USER_FLOATS; ii++ )
    {
        cout << user_float[ii] << " ";
    }
    cout << endl;
    cout << "----------------------------------------------------------" << endl;
}

// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave]
//   0  1  2   3     4  5    6     7   8   9  10
// store a scan with 10 dimensions

GtPlusGadgetImageArray::GtPlusGadgetImageArray() 
:   imageArray_(0)
{
    size_t ii;
    for ( ii=0; ii<GT_DIM_NUM; ii++ )
    {
        matrix_size[ii] = 0;
    }

    max_num_of_images_ = 0;
}

GtPlusGadgetImageArray::GtPlusGadgetImageArray(const GtPlusGadgetImageArray& imArray) : imageArray_(0) 
{
    this->copy(imArray);
}

GtPlusGadgetImageArray::GtPlusGadgetImageArray(size_t aSize[GT_DIM_NUM])
{
    try
    {
        size_t ii;
        for ( ii=0; ii<GT_DIM_NUM; ii++ )
        {
            matrix_size[ii] = aSize[ii];
        }

        size_t len = 1;
        for ( ii=3; ii<GT_DIM_NUM; ii++ )
        {
            len *= matrix_size[ii];
        }

        max_num_of_images_ = len;

        if ( len > 0 )
        {
            imageArray_ = new GtPlusGadgetImageExt[len];
        }
    }
    catch(...)
    {
        GDEBUG_STREAM("Failed in allocate imageArray_" << std::endl);
    }
}

GtPlusGadgetImageArray::~GtPlusGadgetImageArray()
{
    if (imageArray_)
    {
        delete [] imageArray_;
    }
}

void GtPlusGadgetImageArray::resize(size_t aSize[GT_DIM_NUM])
{
    try
    {
        size_t ii;
        for ( ii=0; ii<GT_DIM_NUM; ii++ )
        {
            matrix_size[ii] = aSize[ii];
        }

        size_t len = 1;
        for ( ii=3; ii<GT_DIM_NUM; ii++ )
        {
            len *= matrix_size[ii];
        }

        if ( imageArray_ ) 
        {
            delete [] imageArray_;
            imageArray_ = NULL;
        }

        max_num_of_images_ = len;

        if ( len > 0 )
        {
            imageArray_ = new GtPlusGadgetImageExt[len];
        }
    }
    catch(...)
    {
        GDEBUG_STREAM("Failed in resize GtPlusGadgetImageArray " << std::endl);
    }
}

bool GtPlusGadgetImageArray::copy(const GtPlusGadgetImageArray& imageArray)
{
    try
    {
        if (imageArray_) delete [] imageArray_;
        max_num_of_images_ = 0;

        size_t ii;
        for ( ii=0; ii<GT_DIM_NUM; ii++ )
        {
            matrix_size[ii] = imageArray.matrix_size[ii];
        }

        size_t len = 1;
        for ( ii=3; ii<GT_DIM_NUM; ii++ )
        {
            len *= matrix_size[ii];
        }

        max_num_of_images_ = len;

        if ( len > 0 )
        {
            imageArray_ = new GtPlusGadgetImageExt[len];
        }

        for ( size_t i=0; i<len; i++ )
        {
            imageArray_[i] = imageArray.imageArray_[i];
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusGadgetImageArray::copy(const GtPlusGadgetImageArray& imageArray) ... ");
        return false;
    }

    return true;
}

size_t GtPlusGadgetImageArray::get_offset(size_t slc, size_t e2, size_t con, size_t phs, size_t rep, size_t set, size_t seg, size_t ave)
{
    size_t offset = ave  *matrix_size[9]*matrix_size[8]*matrix_size[7]*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + seg*matrix_size[8]*matrix_size[7]*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + set*matrix_size[7]*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + rep*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + phs*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + con*matrix_size[4]*matrix_size[3]
                    + e2 *matrix_size[3]
                    + slc;
    return offset;
}

// Slice E2 Con Phase Rep Set Seg
void GtPlusGadgetImageArray::findDimIndex(Gadgetron::gtPlus::ISMRMRDDIM& dim, size_t& ind)
{
    switch (dim)
    {
        case Gadgetron::gtPlus::DIM_Slice:
            ind = 3;
        break;

        case Gadgetron::gtPlus::DIM_Encoding2:
            ind = 4;
        break;

        case Gadgetron::gtPlus::DIM_Contrast:
            ind = 5;
        break;

        case Gadgetron::gtPlus::DIM_Phase:
            ind = 6;
        break;

        case Gadgetron::gtPlus::DIM_Repetition:
            ind = 7;
        break;

        case Gadgetron::gtPlus::DIM_Set:
            ind = 8;
        break;

        case Gadgetron::gtPlus::DIM_Segment:
            ind = 9;
        break;

        case Gadgetron::gtPlus::DIM_Average:
            ind = 10;
        break;

        default:
            ind = 0;
    }

    return;
}

bool GtPlusGadgetImageArray::
extractGadgetImageArrayEqual(Gadgetron::gtPlus::ISMRMRDDIM& dim, size_t value, GtPlusGadgetImageArray& imageArray)
{
    try
    {
        size_t dimInd;
        findDimIndex(dim, dimInd);

        GADGET_DEBUG_CHECK_RETURN_FALSE( value >= matrix_size[dimInd] );

        size_t startInd[GT_DIM_NUM-3];
        size_t endInd[GT_DIM_NUM-3];

        for ( size_t d=Gadgetron::gtPlus::DIM_Slice; d<=Gadgetron::gtPlus::DIM_Average; d++ )
        {
            if ( d == dim )
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = value;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = value+1;
            }
            else
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = 0;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = matrix_size[d-Gadgetron::gtPlus::DIM_Slice+3];
            }
        }

        GADGET_CHECK_RETURN_FALSE(getSubImageArray(startInd, endInd, imageArray));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusGadgetImageArray::extractGadgetImageArrayEqual(dim, value) ... ");
        return false;
    }

    return true;
}

bool GtPlusGadgetImageArray::
extractGadgetImageArrayEqual(Gadgetron::gtPlus::ISMRMRDDIM& dim1, size_t value1, Gadgetron::gtPlus::ISMRMRDDIM& dim2, size_t value2, GtPlusGadgetImageArray& imageArray)
{
    try
    {
        size_t dimInd1;
        findDimIndex(dim1, dimInd1);
        GADGET_DEBUG_CHECK_RETURN_FALSE( value1 >= matrix_size[dimInd1] );


        size_t dimInd2;
        findDimIndex(dim2, dimInd2);
        GADGET_DEBUG_CHECK_RETURN_FALSE( value2 >= matrix_size[dimInd2] );

        size_t startInd[GT_DIM_NUM-3];
        size_t endInd[GT_DIM_NUM-3];

        for ( size_t d=Gadgetron::gtPlus::DIM_Slice; d<=Gadgetron::gtPlus::DIM_Average; d++ )
        {
            if ( d == dim1 )
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = value1;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = value1+1;
            }
            else if ( d == dim2 )
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = value2;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = value2+1;
            }
            else
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = 0;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = matrix_size[d-Gadgetron::gtPlus::DIM_Slice+3];
            }
        }

        GADGET_CHECK_RETURN_FALSE(getSubImageArray(startInd, endInd, imageArray));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusGadgetImageArray::extractGadgetImageArrayEqual(dim1, value1, dim2, value2) ... ");
        return false;
    }

    return true;
}

bool GtPlusGadgetImageArray::
extractGadgetImageArrayLessEqual(Gadgetron::gtPlus::ISMRMRDDIM& dim, size_t value, GtPlusGadgetImageArray& imageArray)
{
    try
    {
        size_t dimInd;
        findDimIndex(dim, dimInd);
        GADGET_DEBUG_CHECK_RETURN_FALSE( value >= matrix_size[dimInd] );

        size_t startInd[GT_DIM_NUM-3];
        size_t endInd[GT_DIM_NUM-3];

        for ( size_t d=Gadgetron::gtPlus::DIM_Slice; d<=Gadgetron::gtPlus::DIM_Average; d++ )
        {
            if ( d == dim )
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = 0;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = value+1;
            }
            else
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = 0;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = matrix_size[d-Gadgetron::gtPlus::DIM_Slice+3];
            }
        }

        GADGET_CHECK_RETURN_FALSE(getSubImageArray(startInd, endInd, imageArray));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusGadgetImageArray::extractGadgetImageArrayLessEqual(dim, value) ... ");
        return false;
    }

    return true;
}

bool GtPlusGadgetImageArray::
extractGadgetImageArray_Dim1LessEqual_Dim2Equal(Gadgetron::gtPlus::ISMRMRDDIM& dim1, size_t value1, 
        Gadgetron::gtPlus::ISMRMRDDIM& dim2, size_t value2, GtPlusGadgetImageArray& imageArray)
{
    try
    {
        size_t dimInd1;
        findDimIndex(dim1, dimInd1);

        size_t dimInd2;
        findDimIndex(dim2, dimInd2);

        GADGET_DEBUG_CHECK_RETURN_FALSE( value1 >= matrix_size[dimInd1] );
        GADGET_DEBUG_CHECK_RETURN_FALSE( value2 >= matrix_size[dimInd2] );

        size_t startInd[GT_DIM_NUM];
        size_t endInd[GT_DIM_NUM];

        for ( size_t d=Gadgetron::gtPlus::DIM_Slice; d<=Gadgetron::gtPlus::DIM_Average; d++ )
        {
            if ( d == dim1 )
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = 0;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = value1+1;
            }
            else if ( d == dim2 )
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = value2;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = value2+1;
            }
            else
            {
                startInd[d-Gadgetron::gtPlus::DIM_Slice] = 0;
                endInd[d-Gadgetron::gtPlus::DIM_Slice] = matrix_size[d-Gadgetron::gtPlus::DIM_Slice+3];
            }
        }

        GADGET_CHECK_RETURN_FALSE(getSubImageArray(startInd, endInd, imageArray));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusGadgetImageArray::extractGadgetImageArray_Dim1LessEqual_Dim2Equal(...) ... ");
        return false;
    }

    return true;
}

bool GtPlusGadgetImageArray::
getSubImageArray(size_t* startInd, size_t* endInd, GtPlusGadgetImageArray& imageArray)
{
    try
    {
        size_t aSize[GT_DIM_NUM];
        aSize[0] = matrix_size[0];
        aSize[1] = matrix_size[1];
        aSize[2] = matrix_size[2];

        size_t ii;
        for ( ii=3; ii<GT_DIM_NUM; ii++ )
        {
            aSize[ii] = endInd[ii-3]-startInd[ii-3];
        }

        imageArray.resize(aSize);

        size_t slc, e2, con, phs, rep, set, seg, ave;

        for ( ave=startInd[7]; ave<endInd[7]; ave++ )
        {
            for ( seg=startInd[6]; seg<endInd[6]; seg++ )
            {
                for ( set=startInd[5]; set<endInd[5]; set++ )
                {
                    for ( rep=startInd[4]; rep<endInd[4]; rep++ )
                    {
                        for ( phs=startInd[3]; phs<endInd[3]; phs++ )
                        {
                            for ( con=startInd[2]; con<endInd[2]; con++ )
                            {
                                for ( e2=startInd[1]; e2<endInd[1]; e2++ )
                                {
                                    for ( slc=startInd[0]; slc<endInd[0]; slc++ )
                                    {
                                        size_t offset = this->get_offset(slc, e2, con, phs, rep, set, seg, ave);
                                        size_t offsetDst= imageArray.get_offset(slc-startInd[0], e2-startInd[1], con-startInd[2], phs-startInd[3], rep-startInd[4], set-startInd[5], seg-startInd[6], ave-startInd[7]);

                                        imageArray.imageArray_[offsetDst] = imageArray_[offset];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusGadgetImageArray::getSubImageArray(...) ... ");
        return false;
    }

    return true;
}

void GtPlusGadgetImageArray::dump()
{
    size_t ii;
    GDEBUG_STREAM("GtPlusGadgetImageArray" << std::endl);
    GDEBUG_STREAM("==========================================================" << std::endl);
    GDEBUG_STREAM("matrix_size           : ");
    for ( ii=0; ii<GT_DIM_NUM; ii++ )
    {
        GDEBUG_STREAM(matrix_size[ii] << " ");
    }
    GDEBUG_STREAM(std::endl);
    GDEBUG_STREAM("----------------------------------------------------------" << std::endl);
    if ( imageArray_ )
    {
        int slc, e2, con, phs, rep, set, seg, ave;

        for ( ave=0; ave<matrix_size[10]; ave++ )
        {
            for ( seg=0; seg<matrix_size[9]; seg++ )
            {
                for ( set=0; set<matrix_size[8]; set++ )
                {
                    for ( rep=0; rep<matrix_size[7]; rep++ )
                    {
                        for ( phs=0; phs<matrix_size[6]; phs++ )
                        {
                            for ( con=0; con<matrix_size[5]; con++ )
                            {
                                for ( e2=0; e2<matrix_size[4]; e2++ )
                                {
                                    for ( slc=0; slc<matrix_size[3]; slc++ )
                                    {
                                        size_t offset = get_offset(slc, e2, con, phs, rep, set, seg, ave);
                                        std::cout << "[Slice E2 Contrast Phase Rep Set Seg Ave] = [" 
                                                    << " " << slc 
                                                    << " " << e2 
                                                    << " " << con 
                                                    << " " << phs 
                                                    << " " << rep 
                                                    << " " << set 
                                                    << " " << seg 
                                                    << " " << ave << "]" << std::endl;

                                        imageArray_[offset].dump();
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    GDEBUG_STREAM("==========================================================" << std::endl);
}

}
