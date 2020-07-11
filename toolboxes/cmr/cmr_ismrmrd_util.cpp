/** \file       cmr_ismrmrd_util.cpp
    \brief      Some util functions for ismrmrd.
    \author     Hui Xue
*/

#include "cmr_ismrmrd_util.h"
#include "mri_core_def.h"

namespace Gadgetron
{
    void set_attrib_from_ismrmrd_header(const ISMRMRD::ISMRMRD_ImageHeader& header, ISMRMRD::MetaContainer& attrib)
    {
        attrib.set("PatientPosition", (double)header.position[0]);
        attrib.append("PatientPosition", (double)header.position[1]);
        attrib.append("PatientPosition", (double)header.position[2]);

        attrib.set("read_dir", (double)header.read_dir[0]);
        attrib.append("read_dir", (double)header.read_dir[1]);
        attrib.append("read_dir", (double)header.read_dir[2]);

        attrib.set("phase_dir", (double)header.phase_dir[0]);
        attrib.append("phase_dir", (double)header.phase_dir[1]);
        attrib.append("phase_dir", (double)header.phase_dir[2]);

        attrib.set("slice_dir", (double)header.slice_dir[0]);
        attrib.append("slice_dir", (double)header.slice_dir[1]);
        attrib.append("slice_dir", (double)header.slice_dir[2]);

        attrib.set("patient_table_position", (double)header.patient_table_position[0]);
        attrib.append("patient_table_position", (double)header.patient_table_position[1]);
        attrib.append("patient_table_position", (double)header.patient_table_position[2]);

        attrib.set("FOV", (double)header.field_of_view[0]);
        attrib.append("FOV", (double)header.field_of_view[1]);
        attrib.append("FOV", (double)header.field_of_view[2]);
    }

    bool getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<long>& v)
    {
        try
        {
            size_t num = attrib.length(name.c_str());
            if (num == 0)
            {
                v.clear();
                // GWARN_STREAM("getISMRMRMetaValues, can not find field : " << name);
                return true;
            }

            v.resize(num);

            size_t ii;
            for (ii = 0; ii < num; ii++)
            {
                v[ii] = attrib.as_long(name.c_str(), ii);
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<long>& v) ... ");
            return false;
        }

        return true;
    }

    bool getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<double>& v)
    {
        try
        {
            size_t num = attrib.length(name.c_str());
            if (num == 0)
            {
                v.clear();
                // GWARN_STREAM("getISMRMRMetaValues, can not find field : " << name);
                return true;
            }

            v.resize(num);

            size_t ii;
            for (ii = 0; ii < num; ii++)
            {
                v[ii] = attrib.as_double(name.c_str(), ii);
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<double>& v) ... ");
            return false;
        }

        return true;
    }

    bool getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<std::string>& v)
    {
        try
        {
            size_t num = attrib.length(name.c_str());
            if (num == 0)
            {
                v.clear();
                // GWARN_STREAM("getISMRMRMetaValues, can not find field : " << name);
                return true;
            }

            v.resize(num);

            size_t ii;
            for (ii = 0; ii < num; ii++)
            {
                v[ii] = std::string(attrib.as_str(name.c_str(), ii));
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<std::string>& v) ... ");
            return false;
        }

        return true;
    }

    template <typename T>
    bool setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v)
    {
        try
        {
            size_t num = v.size();
            if (num == 0)
            {
                // GWARN_STREAM("setISMRMRMetaValues, input vector is empty ... " << name);
                return true;
            }

            attrib.set(name.c_str(), v[0]);

            size_t ii;
            for (ii = 1; ii < v.size(); ii++)
            {
                attrib.append(name.c_str(), v[ii]);
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCMR bool setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<long>& v);
    template EXPORTCMR bool setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<double>& v);

    bool setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v)
    {
        try
        {
            size_t num = v.size();
            if (num == 0)
            {
                GWARN_STREAM("setISMRMRMetaValues, input vector is empty ... " << name);
                return true;
            }

            attrib.set(name.c_str(), v[0].c_str());

            size_t ii;
            for (ii = 1; ii < v.size(); ii++)
            {
                attrib.append(name.c_str(), v[ii].c_str());
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v) ... ");
            return false;
        }

        return true;
    }

    template <typename T>
    bool appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v)
    {
        try
        {
            size_t num = v.size();
            if (num == 0)
            {
                GWARN_STREAM("appendISMRMRMetaValues, input vector is empty ... " << name);
                return true;
            }

            attrib.append(name.c_str(), v[0]);

            size_t ii;
            for (ii = 1; ii < v.size(); ii++)
            {
                attrib.append(name.c_str(), v[ii]);
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCMR bool appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<long>& v);
    template EXPORTCMR bool appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<double>& v);

    bool appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v)
    {
        try
        {
            size_t num = v.size();
            if (num == 0)
            {
                GWARN_STREAM("appendISMRMRMetaValues, input vector is empty ... " << name);
                return true;
            }

            attrib.append(name.c_str(), v[0].c_str());

            size_t ii;
            for (ii = 1; ii < v.size(); ii++)
            {
                attrib.append(name.c_str(), v[ii].c_str());
            }
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v) ... ");
            return false;
        }

        return true;
    }

    void create_image_name_from_header(const ISMRMRD::ImageHeader& head, std::string& name)
    {
        try
        {
            std::stringstream os;
            os << "_Series_" << head.image_series_index << "_Image_" << head.image_index
                << "_RO_" << head.matrix_size[0] << "_E1_" << head.matrix_size[1]
                << "_CON_" << head.contrast << "_PHS_" << head.phase << "_REP_" << head.repetition
                << "_SET_" << head.set << "_SLC_" << head.slice << "_AVE_" << head.average;

            name = os.str();
        }
        catch (...)
        {
            GADGET_THROW("Error happened in create_image_name_from_header ... ");
        }
    }

    void create_image_name_from_header(const ISMRMRD::ImageHeader& head, ISMRMRD::MetaContainer& attrib, std::string& name)
    {
        try
        {
            Gadgetron::create_image_name_from_header(head, name);

            if (attrib.length(GADGETRON_DATA_ROLE) > 0)
            {
                std::vector<std::string> data_roles;
                Gadgetron::getISMRMRMetaValues(attrib, GADGETRON_DATA_ROLE, data_roles);

                std::string fname;

                for (size_t ii = 0; ii < data_roles.size(); ii++)
                {
                    fname.append(data_roles[ii]);
                    fname.append("_");
                }

                fname.append(name);
                name = fname;
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in create_image_name_from_header ... ");
        }
    }

    void initialize_image_header_from_acq_header(const ISMRMRD::ISMRMRD_AcquisitionHeader& acq_h, ISMRMRD::ISMRMRD_ImageHeader& im_h)
    {
        size_t ii;

        im_h.version = acq_h.version;
        im_h.data_type = ISMRMRD::ISMRMRD_FLOAT;
        im_h.flags = acq_h.flags;
        im_h.measurement_uid = acq_h.measurement_uid;
        im_h.matrix_size[0] = 1; im_h.matrix_size[1] = 1; im_h.matrix_size[2] = 1;
        im_h.field_of_view[0] = 0; im_h.field_of_view[1] = 0; im_h.field_of_view[2] = 0;
        im_h.channels = 0;

        im_h.position[0] = acq_h.position[0];
        im_h.position[1] = acq_h.position[1];
        im_h.position[2] = acq_h.position[2];

        im_h.read_dir[0] = acq_h.read_dir[0];
        im_h.read_dir[1] = acq_h.read_dir[1];
        im_h.read_dir[2] = acq_h.read_dir[2];

        im_h.phase_dir[0] = acq_h.phase_dir[0];
        im_h.phase_dir[1] = acq_h.phase_dir[1];
        im_h.phase_dir[2] = acq_h.phase_dir[2];

        im_h.slice_dir[0] = acq_h.slice_dir[0];
        im_h.slice_dir[1] = acq_h.slice_dir[1];
        im_h.slice_dir[2] = acq_h.slice_dir[2];

        im_h.patient_table_position[0] = acq_h.patient_table_position[0];
        im_h.patient_table_position[1] = acq_h.patient_table_position[1];
        im_h.patient_table_position[2] = acq_h.patient_table_position[2];

        im_h.average = acq_h.idx.average;
        im_h.slice = acq_h.idx.slice;
        im_h.contrast = acq_h.idx.contrast;
        im_h.phase = acq_h.idx.phase;
        im_h.repetition = acq_h.idx.repetition;
        im_h.set = acq_h.idx.set;
        im_h.acquisition_time_stamp = acq_h.acquisition_time_stamp;

        for (ii = 0; ii < ISMRMRD::ISMRMRD_PHYS_STAMPS; ii++)
            im_h.physiology_time_stamp[ii] = acq_h.physiology_time_stamp[ii];

        im_h.image_type = ISMRMRD::ISMRMRD_IMTYPE_COMPLEX;
        im_h.image_index = 0;
        im_h.image_series_index = 0;

        for (ii = 0; ii < ISMRMRD::ISMRMRD_USER_INTS; ii++)
            im_h.user_int[ii] = acq_h.user_int[ii];

        for (ii = 0; ii < ISMRMRD::ISMRMRD_USER_FLOATS; ii++)
            im_h.user_float[ii] = acq_h.user_float[ii];

        im_h.attribute_string_len = 0;
    }
}
