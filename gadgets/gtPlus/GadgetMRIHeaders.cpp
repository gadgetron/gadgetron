
#include "GadgetMRIHeaders.h"

// --------------------------------------------------------------------

LoopCounters::LoopCounters() 
{
    line = 0;
    acquisition = 0;
    slice = 0;
    partition = 0;
    echo = 0;
    phase = 0;
    repetition = 0;
    set = 0;
    segment = 0;
    channel = 0;
}

LoopCounters::~LoopCounters() {}

void LoopCounters::dump()
{
    std::cout << "[Line Cha Slice Partition Echo Phase Rep Set Seg] = [" 
                    << line 
                    << " " << channel 
                    << " " << slice 
                    << " " << partition 
                    << " " << echo 
                    << " " << phase 
                    << " " << repetition 
                    << " " << set 
                    << " " << segment << "]" << std::endl;
}

// --------------------------------------------------------------------

GadgetMessageAcquisition::GadgetMessageAcquisition() 
{
    flags = 0;
    meas_uid = 0;
    scan_counter = 0;
    time_stamp = 0;
    pmu_time_stamp = 0;
    samples = 0;
    channels = 0;
    centre_column = 0;
    position[0] = 0.0f; position[1] = 0.0f; position[2] = 0.0f;
    quarternion[0] = 1.0f; quarternion[1] = 0.0f; quarternion[2] = 0.0f; quarternion[3] = 0.0f;
    table_position = 0.0f;
}

GadgetMessageAcquisition::~GadgetMessageAcquisition() {}

float GadgetMessageAcquisition::get_position(unsigned int index) 
{
    if (index < 3) 
    {
        return position[index];
    }
    else
    {
        return 0.0f;
    }
}

void GadgetMessageAcquisition::set_position(unsigned int index, float pos)
{
    if (index < 3)
    {
        position[index] = pos;
    }
}

float GadgetMessageAcquisition::get_quarternion(unsigned int index) 
{
    if (index < 4) 
    {
        return quarternion[index];
    }
    else
    {
        return 0.0f;
    }
}

void GadgetMessageAcquisition::set_quarternion(unsigned int index, float quar)
{
    if (index < 4) 
    {
        quarternion[index] = quar;
    }
}

void GadgetMessageAcquisition::dump()
{
    std::cout << "GadgetMessageAcquisition" << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    std::cout << "flags            : " << flags << std::endl;
    std::cout << "meas_uid         : " << meas_uid << std::endl;
    std::cout << "scan_counter     : " << scan_counter << std::endl;
    std::cout << "time_stamp       : " << time_stamp << std::endl;
    std::cout << "pmu_time_stamp   : " << pmu_time_stamp << std::endl;
    std::cout << "samples          : " << samples << std::endl;
    std::cout << "channels         : " << channels << std::endl;
    std::cout << "position         : " << position[0] << " " << position[1] << " " << position[2] << std::endl;
    std::cout << "quarternion      : " << quarternion[0] << " " << quarternion[1] << " " << quarternion[2] << " " << quarternion[3] << std::endl;
    std::cout << "table_position   : " << table_position << std::endl;
    std::cout << "idx     : ";            idx.dump();
    std::cout << "min_idx : ";            min_idx.dump();
    std::cout << "max_idx : ";            max_idx.dump();
    std::cout << "----------------------------------------------------------" << std::endl;
}

// --------------------------------------------------------------------

GadgetMessageImage::GadgetMessageImage()
{
    flags = 0;

    matrix_size[0] = 0;
    matrix_size[1] = 0;
    matrix_size[2] = 0;

    channels = 0;

    position[0] = 0.0f;
    position[1] = 0.0f;
    position[2] = 0.0f;

    quarternion[0] = 1.0f;
    quarternion[1] = 0.0f;
    quarternion[2] = 0.0f;
    quarternion[3] = 0.0f;

    table_position = 0.0f;

    time_stamp = 0;
    pmu_time_stamp = 0;
    image_format = 0;
    image_type = 0;
    image_index = 0;
    image_series_index = 0;
}

GadgetMessageImage::~GadgetMessageImage() {}

void GadgetMessageImage::copy(GadgetMessageImage& aMessageImage)
{
    flags = aMessageImage.flags;

    matrix_size[0] = aMessageImage.matrix_size[0];
    matrix_size[1] = aMessageImage.matrix_size[1];
    matrix_size[2] = aMessageImage.matrix_size[2];

    channels = aMessageImage.channels;

    position[0] = aMessageImage.position[0];
    position[1] = aMessageImage.position[1];
    position[2] = aMessageImage.position[2];

    quarternion[0] = aMessageImage.quarternion[0];
    quarternion[1] = aMessageImage.quarternion[1];
    quarternion[2] = aMessageImage.quarternion[2];
    quarternion[3] = aMessageImage.quarternion[3];

    table_position = aMessageImage.table_position;

    time_stamp = aMessageImage.time_stamp;
    pmu_time_stamp = aMessageImage.pmu_time_stamp;
    image_format = aMessageImage.image_format;
    image_type = aMessageImage.image_type;
    image_index = aMessageImage.image_index;
    image_series_index = aMessageImage.image_series_index;
}

ACE_UINT16 GadgetMessageImage::get_matrix_size(unsigned int index) 
{
    if (index < 3) 
    {
        return matrix_size[index];
    }
    else
    {
        return 0;
    }
}

void GadgetMessageImage::set_matrix_size(unsigned int index, ACE_UINT16 size)
{
    if (index < 3) 
    {
        matrix_size[index] = size;
    }
}

float GadgetMessageImage::get_position(unsigned int index) 
{
    if (index < 3) 
    {
        return position[index];
    }
    else
    {
        return 0.0f;
    }
}

void GadgetMessageImage::set_position(unsigned int index, float pos)
{
    if (index < 3)
    {
        position[index] = pos;
    }
}

float GadgetMessageImage::get_quarternion(unsigned int index)
{
    if (index < 4)
    {
        return quarternion[index];
    }
    else
    {
        return 0.0f;
    }
}

void GadgetMessageImage::set_quarternion(unsigned int index, float quar)
{
    if (index < 4)
    {
        quarternion[index] = quar;
    }
}

void GadgetMessageImage::dumpInfo()
{
    std::cout << "flags                 : " << flags << std::endl;
    std::cout << "matrix_size           : " << matrix_size[0] << " " << matrix_size[1] << " " << matrix_size[2] << std::endl;
    std::cout << "channels              : " << channels << std::endl;
    std::cout << "position              : " << position[0] << " " << position[1] << " " << position[2] << std::endl;
    std::cout << "quarternion           : " << quarternion[0] << " " << quarternion[1] << " " << quarternion[2] << " " << quarternion[3] << std::endl;
    std::cout << "table_position        : " << table_position << std::endl;
    std::cout << "data_idx_min          : ";   data_idx_min.dump();
    std::cout << "data_idx_max          : ";   data_idx_max.dump();
    std::cout << "data_idx_current      : ";   data_idx_current.dump();
    std::cout << "time_stamp            : " << time_stamp << std::endl;
    std::cout << "pmu_time_stamp        : " << pmu_time_stamp << std::endl;
    std::cout << "image_format          : " << image_format << std::endl;
    std::cout << "image_type            : " << image_type << std::endl;
    std::cout << "image_index           : " << image_index << std::endl;
    std::cout << "image_series_index    : " << image_series_index << std::endl;
}

void GadgetMessageImage::dump()
{
    std::cout << "GadgetMessageImage" << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    dumpInfo();
    std::cout << "----------------------------------------------------------" << std::endl;
}
