
#include "GadgetMRIHeadersExt.h"
#include "GadgetIsmrmrdReadWrite.h"
// #include <iostream>

// --------------------------------------------------------------------

GadgetMessageImageExt::GadgetMessageImageExt() : ISMRMRD::ImageHeader()
{
    time_stamps.clear();
    pmu_time_stamps.clear();
}

GadgetMessageImageExt::~GadgetMessageImageExt() { }

void GadgetMessageImageExt::set_matrix_size(unsigned int index, ACE_UINT16 size)
{
    if (index < 3) 
    {
        matrix_size[index] = size;
    }

    if ( index == 1 )
    {
        time_stamps.clear();
        time_stamps.resize(matrix_size[1], -1);
        pmu_time_stamps.clear();
        pmu_time_stamps.resize(matrix_size[1], -1);
    }
}

void GadgetMessageImageExt::copy(GadgetMessageImageExt& aMessageImage)
{
    flags = aMessageImage.flags;

    matrix_size[0] = aMessageImage.matrix_size[0];
    matrix_size[1] = aMessageImage.matrix_size[1];
    matrix_size[2] = aMessageImage.matrix_size[2];

    channels = aMessageImage.channels;

    position[0] = aMessageImage.position[0];
    position[1] = aMessageImage.position[1];
    position[2] = aMessageImage.position[2];

    read_dir[0] = aMessageImage.read_dir[0];
    read_dir[1] = aMessageImage.read_dir[1];
    read_dir[2] = aMessageImage.read_dir[2];

    phase_dir[0] = aMessageImage.phase_dir[0];
    phase_dir[1] = aMessageImage.phase_dir[1];
    phase_dir[2] = aMessageImage.phase_dir[2];

    slice_dir[0] = aMessageImage.slice_dir[0];
    slice_dir[1] = aMessageImage.slice_dir[1];
    slice_dir[2] = aMessageImage.slice_dir[2];

    patient_table_position[0] = aMessageImage.patient_table_position[0];
    patient_table_position[1] = aMessageImage.patient_table_position[1];
    patient_table_position[2] = aMessageImage.patient_table_position[2];

    acquisition_time_stamp = aMessageImage.acquisition_time_stamp;

    physiology_time_stamp[0] = aMessageImage.physiology_time_stamp[0];
    physiology_time_stamp[1] = aMessageImage.physiology_time_stamp[1];
    physiology_time_stamp[2] = aMessageImage.physiology_time_stamp[2];

    image_data_type = aMessageImage.image_data_type;
    image_type = aMessageImage.image_type;
    image_index = aMessageImage.image_index;
    image_series_index = aMessageImage.image_series_index;

    memcpy(user_int, aMessageImage.user_int, sizeof(int32_t)*ISMRMRD_USER_INTS);
    memcpy(user_float, aMessageImage.user_float, sizeof(float)*ISMRMRD_USER_FLOATS);

    time_stamps = aMessageImage.time_stamps;
    pmu_time_stamps = aMessageImage.pmu_time_stamps;
}

void GadgetMessageImageExt::dump()
{
    GDEBUG_STREAM("GadgetMessageImageExt" << std::endl);
    GDEBUG_STREAM("----------------------------------------------------------" << std::endl);
    //dumpInfo();
    GDEBUG_STREAM("----------------------------------------------------------" << std::endl);
}

// --------------------------------------------------------------------

// [Col Line Cha Slice Partition Echo Phase Rep Set Seg]
//   0   1    2   3     4         5    6     7   8   9
// store a scan with 10 dimensions
GadgetMessageImageArray::GadgetMessageImageArray() 
:   imageArray_(0),
    kSpace_centre_col_no(0), 
    kSpace_centre_line_no(0), 
    kSpace_centre_partition_no(0), 
    kSpace_max_acquired_col_no(0), 
    kSpace_max_acquired_line_no(0), 
    kSpace_max_acquired_partition_no(0)
{

}

GadgetMessageImageArray::GadgetMessageImageArray(int aSize[10])
{
    try
    {
        unsigned int ii;
        for ( ii=0; ii<10; ii++ )
        {
            matrix_size[ii] = aSize[ii];
        }

        unsigned int len = 1;
        for ( ii=3; ii<10; ii++ )
        {
            len *= matrix_size[ii];
        }

        if ( len > 0 )
        {
            imageArray_ = new GadgetMessageImageExt[len];
        }

        kSpace_centre_col_no = matrix_size[0]/2;
        kSpace_centre_line_no = matrix_size[1]/2;
        kSpace_centre_partition_no = matrix_size[4]/2;

        kSpace_max_acquired_col_no = matrix_size[0]-1;
        kSpace_max_acquired_line_no = matrix_size[1]-1;
        kSpace_max_acquired_partition_no = matrix_size[4]-1;
    }
    catch(...)
    {
        GDEBUG_STREAM("Failed in allocate imageArray_" << std::endl);
    }
}

GadgetMessageImageArray::~GadgetMessageImageArray()
{
    if (imageArray_)
    {
        delete [] imageArray_;
    }
}

void GadgetMessageImageArray::resize(int aSize[10])
{
    try
    {
        unsigned int ii;
        for ( ii=0; ii<10; ii++ )
        {
            matrix_size[ii] = aSize[ii];
        }

        unsigned int len = 1;
        for ( ii=3; ii<10; ii++ )
        {
            len *= matrix_size[ii];
        }

        if ( imageArray_ ) 
        {
            delete [] imageArray_;
            imageArray_ = NULL;
        }

        if ( len > 0 )
        {
            imageArray_ = new GadgetMessageImageExt[len];
        }

        kSpace_centre_col_no = matrix_size[0]/2;
        kSpace_centre_line_no = matrix_size[1]/2;
        kSpace_centre_partition_no = matrix_size[4]/2;

        kSpace_max_acquired_col_no = matrix_size[0]-1;
        kSpace_max_acquired_line_no = matrix_size[1]-1;
        kSpace_max_acquired_partition_no = matrix_size[4]-1;
    }
    catch(...)
    {
        GDEBUG_STREAM("Failed in resize GadgetMessageImageArray " << std::endl);
    }
}

void GadgetMessageImageArray::copy(GadgetMessageImageArray& imageArray)
{
    if (imageArray_) delete [] imageArray_;

    unsigned int ii;
    for ( ii=0; ii<10; ii++ )
    {
        matrix_size[ii] = imageArray.matrix_size[ii];
    }

    unsigned int len = 1;
    for ( ii=3; ii<10; ii++ )
    {
        len *= matrix_size[ii];
    }

    kSpace_centre_col_no = imageArray.kSpace_centre_col_no;
    kSpace_centre_line_no = imageArray.kSpace_centre_line_no;
    kSpace_centre_partition_no = imageArray.kSpace_centre_partition_no;

    kSpace_max_acquired_col_no = imageArray.kSpace_max_acquired_col_no;
    kSpace_max_acquired_line_no = imageArray.kSpace_max_acquired_line_no;
    kSpace_max_acquired_partition_no = imageArray.kSpace_max_acquired_partition_no;

    if ( len > 0 )
    {
        imageArray_ = new GadgetMessageImageExt[len];
    }

    for ( unsigned int i=0; i<len; i++ )
    {
        imageArray_[i] = imageArray.imageArray_[i];
    }
}

int GadgetMessageImageArray::get_offset(int slc, int par, int eco, int phs, int rep, int set, int seg)
{
    int offset = seg*matrix_size[8]*matrix_size[7]*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + set*matrix_size[7]*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + rep*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + phs*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + eco*matrix_size[4]*matrix_size[3]
                    + par*matrix_size[3]
                    + slc;
    return offset;
}

void GadgetMessageImageArray::extractMessageImageArrayForSLC(int slc, GadgetMessageImageArray& imageArray)
{
    if ( slc >= matrix_size[3] )
    {
        GDEBUG_STREAM("extractMessageImageArrayForSLC error - slc >= matrix_size[3] " << std::endl);
        return;
    }

    int aSize[10];

    unsigned int ii;
    for ( ii=0; ii<10; ii++ )
    {
        aSize[ii] = matrix_size[ii];
    }

    aSize[3] = 1;

    imageArray.resize(aSize);

    imageArray.kSpace_centre_col_no = kSpace_centre_col_no;
    imageArray.kSpace_centre_line_no = kSpace_centre_line_no;
    imageArray.kSpace_centre_partition_no = kSpace_centre_partition_no;
    imageArray.kSpace_max_acquired_col_no = kSpace_max_acquired_col_no;
    imageArray.kSpace_max_acquired_line_no = kSpace_max_acquired_line_no;
    imageArray.kSpace_max_acquired_partition_no = kSpace_max_acquired_partition_no;

    int par, eco, phs, rep, set, seg;

    int PAR = matrix_size[4];
    int ECO = matrix_size[5];
    int PHS = matrix_size[6];
    int REP = matrix_size[7];
    int SET = matrix_size[8];
    int SEG = matrix_size[9];

    for ( seg=0; seg<SEG; seg++ )
    {
        for ( set=0; set<SET; set++ )
        {
            for ( rep=0; rep<REP; rep++ )
            {
                for ( phs=0; phs<PHS; phs++ )
                {
                    for ( eco=0; eco<ECO; eco++ )
                    {
                        for ( par=0; par<PAR; par++ )
                        {
                            int offset = this->get_offset(slc, par, eco, phs, rep, set, seg);
                            int offsetSLC = imageArray.get_offset(0, par, eco, phs, rep, set, seg);

                            imageArray.imageArray_[offsetSLC] = imageArray_[offset];
                        }
                    }
                }
            }
        }
    }
}

void GadgetMessageImageArray::extractMessageImageArrayForREP(int rep, GadgetMessageImageArray& imageArray)
{
    if ( rep >= matrix_size[7] )
    {
        GDEBUG_STREAM("extractMessageImageArrayForSLC error - rep >= matrix_size[7] " << std::endl);
        return;
    }

    int aSize[10];

    unsigned int ii;
    for ( ii=0; ii<10; ii++ )
    {
        aSize[ii] = matrix_size[ii];
    }

    aSize[7] = 1;

    imageArray.resize(aSize);

    imageArray.kSpace_centre_col_no = kSpace_centre_col_no;
    imageArray.kSpace_centre_line_no = kSpace_centre_line_no;
    imageArray.kSpace_centre_partition_no = kSpace_centre_partition_no;
    imageArray.kSpace_max_acquired_col_no = kSpace_max_acquired_col_no;
    imageArray.kSpace_max_acquired_line_no = kSpace_max_acquired_line_no;
    imageArray.kSpace_max_acquired_partition_no = kSpace_max_acquired_partition_no;

    int par, eco, phs, slc, set, seg;

    int SLC = matrix_size[3];
    int PAR = matrix_size[4];
    int ECO = matrix_size[5];
    int PHS = matrix_size[6];
    int SET = matrix_size[8];
    int SEG = matrix_size[9];

    for ( seg=0; seg<SEG; seg++ )
    {
        for ( set=0; set<SET; set++ )
        {
            for ( slc=0; slc<SLC; slc++ )
            {
                for ( phs=0; phs<PHS; phs++ )
                {
                    for ( eco=0; eco<ECO; eco++ )
                    {
                        for ( par=0; par<PAR; par++ )
                        {
                            int offset = this->get_offset(slc, par, eco, phs, rep, set, seg);
                            int offsetREP = imageArray.get_offset(slc, par, eco, phs, 0, set, seg);

                            imageArray.imageArray_[offsetREP] = imageArray_[offset];
                        }
                    }
                }
            }
        }
    }
}

void GadgetMessageImageArray::dump()
{
    unsigned int ii;
    GDEBUG_STREAM("GadgetMessageImageArray" << std::endl);
    GDEBUG_STREAM("==========================================================" << std::endl);
    GDEBUG_STREAM("matrix_size           : ");
    for ( ii=0; ii<10; ii++ )
    {
        GDEBUG_STREAM(matrix_size[ii] << " ");
    }
    GDEBUG_STREAM(std::endl);
    GDEBUG_STREAM("----------------------------------------------------------" << std::endl);
    GDEBUG_STREAM("kSpace_centre_col_no             : " << kSpace_centre_col_no << std::endl);
    GDEBUG_STREAM("kSpace_max_acquired_col_no       : " << kSpace_max_acquired_col_no << std::endl);
    GDEBUG_STREAM("----------------------------------------------------------" << std::endl);
    GDEBUG_STREAM("kSpace_centre_line_no            : " << kSpace_centre_line_no << std::endl);
    GDEBUG_STREAM("kSpace_max_acquired_line_no      : " << kSpace_max_acquired_line_no << std::endl);
    GDEBUG_STREAM("----------------------------------------------------------" << std::endl);
    GDEBUG_STREAM("kSpace_centre_partition_no       : " << kSpace_centre_partition_no << std::endl);
    GDEBUG_STREAM("kSpace_max_acquired_partition_no : " << kSpace_max_acquired_partition_no << std::endl);
    GDEBUG_STREAM("----------------------------------------------------------" << std::endl);
    if ( imageArray_ )
    {
        int slc, par, eco, phs, rep, set, seg;
        for ( seg=0; seg<matrix_size[9]; seg++ )
        {
            for ( set=0; set<matrix_size[8]; set++ )
            {
                for ( rep=0; rep<matrix_size[7]; rep++ )
                {
                    for ( phs=0; phs<matrix_size[6]; phs++ )
                    {
                        for ( eco=0; eco<matrix_size[5]; eco++ )
                        {
                            for ( par=0; par<matrix_size[4]; par++ )
                            {
                                for ( slc=0; slc<matrix_size[3]; slc++ )
                                {
                                    int offset = get_offset(slc, par, eco, phs, rep, set, seg);
                                    std::cout << "[Slice Partition Echo Phase Rep Set Seg] = [" 
                                                << " " << slc 
                                                << " " << par 
                                                << " " << eco 
                                                << " " << phs 
                                                << " " << rep 
                                                << " " << set 
                                                << " " << seg << "]" << std::endl;

                                    imageArray_[offset].dump();
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

// --------------------------------------------------------------------

KSpaceBuffer::KSpaceBuffer() 
: isIPAT(false) 
{

}

KSpaceBuffer::~KSpaceBuffer()
{

}
