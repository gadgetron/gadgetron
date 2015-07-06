
#include "GtPlusReconGadgetUtil.h"

#include <boost/filesystem.hpp>
using namespace boost::filesystem;

namespace Gadgetron
{

    bool findCalibMode(ISMRMRD::IsmrmrdHeader& h, Gadgetron::ISMRMRDCALIBMODE& CalibMode, ISMRMRDDIM& InterleaveDim, double& acceFactorE1, double& acceFactorE2, bool verbose)
    {
        try
        {
            if (!h.encoding[0].parallelImaging)
            {
                GERROR_STREAM("Parallel Imaging section not found in header");
                return false;
            }

            ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

            acceFactorE1 = (double)(p_imaging.accelerationFactor.kspace_encoding_step_1);
            acceFactorE2 = (double)(p_imaging.accelerationFactor.kspace_encoding_step_2);

            GDEBUG_CONDITION_STREAM(verbose, "acceFactorE1 is " << acceFactorE1);
            GDEBUG_CONDITION_STREAM(verbose, "acceFactorE2 is " << acceFactorE2);

            if ( !p_imaging.calibrationMode.is_present() )
            {
                GERROR_STREAM("Parallel calibration mode not found in header");
                return false;
            }

            std::string calib = *p_imaging.calibrationMode;
            if ( calib.compare("interleaved") == 0 )
            {
                CalibMode = Gadgetron::ISMRMRD_interleaved;
                GDEBUG_CONDITION_STREAM(verbose, "Calibration mode is interleaved");

                if ( p_imaging.interleavingDimension )
                {
                    if ( p_imaging.interleavingDimension->compare("phase") == 0 )
                    {
                        InterleaveDim = Gadgetron::DIM_Phase;
                    }
                    else if ( p_imaging.interleavingDimension->compare("repetition") == 0 )
                    {
                        InterleaveDim = Gadgetron::DIM_Repetition;
                    }
                    else if ( p_imaging.interleavingDimension->compare("average") == 0 )
                    {
                        InterleaveDim = Gadgetron::DIM_Average;
                    }
                    else if ( p_imaging.interleavingDimension->compare("contrast") == 0 )
                    {
                        InterleaveDim = Gadgetron::DIM_Contrast;
                    }
                    else if ( p_imaging.interleavingDimension->compare("other") == 0 )
                    {
                        InterleaveDim = Gadgetron::DIM_other1;
                    }
                    else
                    {
                        GERROR_STREAM("Unknown interleaving dimension. Bailing out");
                        return false;
                    }
                }
            }
            else if ( calib.compare("embedded") == 0 )
            {
                CalibMode = Gadgetron::ISMRMRD_embedded;
                GDEBUG_CONDITION_STREAM(verbose, "Calibration mode is embedded");
            }
            else if ( calib.compare("separate") == 0 )
            {
                CalibMode = Gadgetron::ISMRMRD_separate;
                GDEBUG_CONDITION_STREAM(verbose, "Calibration mode is separate");
            }
            else if ( calib.compare("external") == 0 )
            {
                CalibMode = Gadgetron::ISMRMRD_external;
            }
            else if ( (calib.compare("other") == 0) && acceFactorE1==1 && acceFactorE2==1 )
            {
                CalibMode = Gadgetron::ISMRMRD_noacceleration;
                acceFactorE1=1;
            }
            else if ( (calib.compare("other") == 0) &&  (acceFactorE1>1 || acceFactorE2>1) )
            {
                CalibMode = Gadgetron::ISMRMRD_interleaved;
                acceFactorE1=2;
                InterleaveDim = Gadgetron::DIM_Phase;
            }
            else
            {
                GERROR_STREAM("Failed to process parallel imaging calibration mode");
                return false;
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in findCalibMode(...) ... ");
            return false;
        }

        return true;
    }

    bool findEncodingLimits(ISMRMRD::IsmrmrdHeader& h, ISMRMRD::EncodingCounters& meas_max_idx, bool verbose)
    {
        try
        {
            ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
            ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
            ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

            meas_max_idx.kspace_encode_step_1 = (uint16_t)e_space.matrixSize.y-1;

            meas_max_idx.set = (e_limits.set && (e_limits.set->maximum>0)) ? e_limits.set->maximum : 0;
            meas_max_idx.phase = (e_limits.phase && (e_limits.phase->maximum>0)) ? e_limits.phase->maximum : 0;

            meas_max_idx.kspace_encode_step_2 = (uint16_t)e_space.matrixSize.z-1;

            meas_max_idx.contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum : 0;

            meas_max_idx.slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;

            meas_max_idx.repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;

            meas_max_idx.average = e_limits.average ? e_limits.average->maximum : 0;

            // always combine the SEG
            meas_max_idx.segment = 0;
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in findEncodingLimits(...) ... ");
            return false;
        }

        return true;
    }

    void findMatrixSizeEncoding(ISMRMRD::IsmrmrdHeader& h, size_t matrix_size_encoding[3])
    {
        matrix_size_encoding[0] = h.encoding[0].encodedSpace.matrixSize.x;
        matrix_size_encoding[1] = h.encoding[0].encodedSpace.matrixSize.y;
        matrix_size_encoding[2] = h.encoding[0].encodedSpace.matrixSize.z;
    }

    void findFOVEncoding(ISMRMRD::IsmrmrdHeader& h, float field_of_view_encoding[3])
    {
        field_of_view_encoding[0] = h.encoding[0].encodedSpace.fieldOfView_mm.x;
        field_of_view_encoding[1] = h.encoding[0].encodedSpace.fieldOfView_mm.y;
        field_of_view_encoding[2] = h.encoding[0].encodedSpace.fieldOfView_mm.z;
    }

    void findMatrixSizeRecon(ISMRMRD::IsmrmrdHeader& h, size_t matrix_size_recon[3])
    {
        matrix_size_recon[0] = h.encoding[0].reconSpace.matrixSize.x;
        matrix_size_recon[1] = h.encoding[0].reconSpace.matrixSize.y;
        matrix_size_recon[2] = h.encoding[0].reconSpace.matrixSize.z;
    }

    void findFOVRecon(ISMRMRD::IsmrmrdHeader& h, float field_of_view_recon[3])
    {
        field_of_view_recon[0] = h.encoding[0].reconSpace.fieldOfView_mm.x;
        field_of_view_recon[1] = h.encoding[0].reconSpace.fieldOfView_mm.y;
        field_of_view_recon[2] = h.encoding[0].reconSpace.fieldOfView_mm.z;
    }

    bool checkReadoutStatus(uint64_t flag, int samples, Gadgetron::ISMRMRDCALIBMODE& CalibMode, int roLen, 
        bool& bIsKSpace, bool& bIsRef, bool& bIsNoise, 
        bool& bIsPhaseCorr, bool& bIsReflect, bool& bIsOther, 
        bool& bIsNavigator, bool& bIsRTFeedback, bool& bIsHPFeedback, 
        bool& bIsDummyScan)
    {
        try
        {
            bIsNoise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(flag);
            bool is_ref = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(flag);
            bool is_ref_kspace = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(flag);
            bIsReflect = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_REVERSE).isSet(flag);
            bIsPhaseCorr = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PHASECORR_DATA).isSet(flag);
            bIsNavigator = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NAVIGATION_DATA).isSet(flag);
            bIsRTFeedback = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_RTFEEDBACK_DATA).isSet(flag);
            bIsHPFeedback = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_HPFEEDBACK_DATA).isSet(flag);
            bIsDummyScan = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_DUMMYSCAN_DATA).isSet(flag);

            bIsKSpace = false;
            bIsRef = false;
            bIsOther = false;

            if ( bIsNoise || bIsDummyScan )
            {
                return true;
            }

            if ( CalibMode==ISMRMRD_noacceleration )
            {
                bIsKSpace = true;
                bIsRef = false;
            }

            // in interleaved mode, only store the image data
            if ( CalibMode==ISMRMRD_interleaved )
            {
                bIsKSpace = true;
                bIsRef = false;
            }

            // in embedded, kspace stores only the undersampled lines
            // ref stores all lines used for references
            if ( CalibMode==ISMRMRD_embedded )
            {
                if ( is_ref && !is_ref_kspace )
                {
                    bIsKSpace = false;
                    bIsRef = true;
                }

                if ( !is_ref && is_ref_kspace )
                {
                    bIsKSpace = true;
                    bIsRef = true;
                }

                if ( is_ref && is_ref_kspace )
                {
                    bIsKSpace = true;
                    bIsRef = true;
                }

                if ( !is_ref && !is_ref_kspace )
                {
                    bIsKSpace = true;
                    bIsRef = false;
                }
            }

            // in separate mode
            if ( CalibMode==ISMRMRD_separate 
                || CalibMode==ISMRMRD_external )
            {
                if ( is_ref )
                {
                    bIsKSpace = false;
                    bIsRef = true;
                }

                if ( !is_ref )
                {
                    bIsKSpace = true;
                    bIsRef = false;
                }
            }

            // store other data, e.g. AIF
            // only for tpat
            if ( !is_ref && !is_ref_kspace && (samples!=roLen) )
            {
                bIsOther = true;
                bIsKSpace = false;
                bIsRef = false;
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in checkReadoutStatus(...) ... ");
            return false;
        }

        return true;
    }

    bool estimateMaxSEGForRetroGating(Gadgetron::ISMRMRDCALIBMODE CalibMode, 
        double acceFactorE1, double acceFactorE2, 
        size_t retro_gated_segment_size, 
        uint16_t E1, uint16_t embedded_ref_lines_E1, 
        uint16_t E2, uint16_t embedded_ref_lines_E2, 
        uint16_t& segment, bool verbose)
    {
        try
        {
            if ( acceFactorE2 <= 1 )
            {
                if ( CalibMode == ISMRMRD_embedded )
                {
                    segment = (uint16_t)std::ceil( (double)E1/acceFactorE1/retro_gated_segment_size 
                        + (acceFactorE1-1)*(double)embedded_ref_lines_E1/acceFactorE1/retro_gated_segment_size );
                }
                else
                {
                    segment = (uint16_t)std::ceil( (double)E1/acceFactorE1/retro_gated_segment_size );
                }
            }
            else
            {
                if ( CalibMode == ISMRMRD_embedded )
                {
                    segment = (uint16_t)std::ceil( (double)E1*E2/(acceFactorE1*acceFactorE2*retro_gated_segment_size) 
                        + (acceFactorE1*acceFactorE2-1)*(double)(embedded_ref_lines_E1*embedded_ref_lines_E2)/(acceFactorE1*acceFactorE2*retro_gated_segment_size) );
                }
                else
                {
                    segment = (uint16_t)std::ceil( (double)E1*E2/(acceFactorE1*acceFactorE2*retro_gated_segment_size) );
                }
            }

            if ( segment > 1 ) segment--;
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in estimateMaxSEGForRetroGating(...) ... ");
            return false;
        }

        return true;
    }

    void getDebugFolderPath(const std::string& debugFolder, std::string& debugFolderPath, bool verbose)
    {
      debugFolderPath = getenv("GADGETRON_DEBUG_FOLDER");
      if ( debugFolderPath.empty() )
      {
#ifdef _WIN32
            debugFolderPath = "c:/temp/gadgetron";
#else
            debugFolderPath = "/tmp/gadgetron";
#endif // _WIN32
        }

        debugFolderPath.append("/");
        debugFolderPath.append(debugFolder);
        debugFolderPath.append("/");

        createFolderWithAllPermissions(debugFolderPath);

        GDEBUG_CONDITION_STREAM(verbose, "Debug folder is " << debugFolderPath);
    }

    bool createFolderWithAllPermissions(const std::string& workingdirectory)
    {
        if ( !boost::filesystem::exists(workingdirectory) )
        {
            boost::filesystem::path workingPath(workingdirectory);
            if ( !boost::filesystem::create_directory(workingPath) )
            {
	      GERROR("Error creating the working directory.\n");
	      return false;
            }

            // set the permission for the folder
#ifdef _WIN32
            try
            {
                boost::filesystem::permissions(workingPath, all_all);
            }
            catch(...)
            {
	      GERROR("Error changing the permission of the working directory.\n");
	      return false;
            }
#else
            // in case an older version of boost is used in non-win system
            // the system call is used
            int res = chmod(workingPath.string().c_str(), S_IRUSR|S_IWUSR|S_IXUSR|S_IRGRP|S_IWGRP|S_IXGRP|S_IROTH|S_IWOTH|S_IXOTH);
            if ( res != 0 )
            {
	      GERROR("Error changing the permission of the working directory.\n");
	      return false;
            }
#endif // _WIN32
        }

        return true;
    }

    bool getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<long>& v)
    {
        try
        {
            size_t num = attrib.length(name.c_str());
            if ( num == 0 )
            {
                v.clear();
                GWARN_STREAM("getISMRMRMetaValues, can not find field : " << name);
                return true;
            }

            v.resize(num);

            size_t ii;
            for ( ii=0; ii<num; ii++ )
            {
                v[ii] = attrib.as_long(name.c_str(), ii);
            }
        }
        catch(...)
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
            if ( num == 0 )
            {
                v.clear();
                GWARN_STREAM("getISMRMRMetaValues, can not find field : " << name);
                return true;
            }

            v.resize(num);

            size_t ii;
            for ( ii=0; ii<num; ii++ )
            {
                v[ii] = attrib.as_double(name.c_str(), ii);
            }
        }
        catch(...)
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
            if ( num == 0 )
            {
                v.clear();
                GWARN_STREAM("getISMRMRMetaValues, can not find field : " << name);
                return true;
            }

            v.resize(num);

            size_t ii;
            for ( ii=0; ii<num; ii++ )
            {
                v[ii] = std::string( attrib.as_str(name.c_str(), ii) );
            }
        }
        catch(...)
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
            if ( num == 0 )
            {
                GWARN_STREAM("setISMRMRMetaValues, input vector is empty ... " << name);
                return true;
            }

            attrib.set(name.c_str(), v[0]);

            size_t ii;
            for ( ii=1; ii<v.size(); ii++ )
            {
                attrib.append(name.c_str(), v[ii]);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v) ... ");
            return false;
        }

        return true;
    }

    template EXPORTGTPLUSGADGET bool setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<long>& v);
    template EXPORTGTPLUSGADGET bool setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<double>& v);

    bool setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v)
    {
        try
        {
            size_t num = v.size();
            if ( num == 0 )
            {
                GWARN_STREAM("setISMRMRMetaValues, input vector is empty ... " << name);
                return true;
            }

            attrib.set(name.c_str(), v[0].c_str());

            size_t ii;
            for ( ii=1; ii<v.size(); ii++ )
            {
                attrib.append(name.c_str(), v[ii].c_str());
            }
        }
        catch(...)
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
            if ( num == 0 )
            {
                GWARN_STREAM("appendISMRMRMetaValues, input vector is empty ... " << name);
                return true;
            }

            attrib.append(name.c_str(), v[0]);

            size_t ii;
            for ( ii=1; ii<v.size(); ii++ )
            {
                attrib.append(name.c_str(), v[ii]);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v) ... ");
            return false;
        }

        return true;
    }

    template EXPORTGTPLUSGADGET bool appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<long>& v);
    template EXPORTGTPLUSGADGET bool appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<double>& v);

    bool appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v)
    {
        try
        {
            size_t num = v.size();
            if ( num == 0 )
            {
                GWARN_STREAM("appendISMRMRMetaValues, input vector is empty ... " << name);
                return true;
            }

            attrib.append(name.c_str(), v[0].c_str());

            size_t ii;
            for ( ii=1; ii<v.size(); ii++ )
            {
                attrib.append(name.c_str(), v[ii].c_str());
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v) ... ");
            return false;
        }

        return true;
    }

    bool PatientCoordinateSystemToDeviceCoordinateSystem(double& x, double& y, double& z, const std::string& position)
    {
        // this is following dicom tag (0020, 0037)

        if ( position == "HFS" ) // Head-first supine (HFS)
        {
            y = -y;
            z = -z;
        }
        else if ( position == "HFP" ) // Head-first prone (HFP)
        {
            x = -x;
            z = -z;
        }
        else if ( position == "HFDR" ) // Head-first decubitus-right 
        {
            double v = x;
            x = y;
            y = v;
            z = -z;
        }
        else if ( position == "HFDL" ) // Head-first decubitus-left (HFDL)
        {
            double v = x;
            x = y;
            y = v;

            x = -x;
            y = -y;
            z = -z;
        }
        else if ( position == "FFDR" ) // Feet-first decubitus-right (FFDR)
        {
            double v = x;
            x = y;
            y = v;

            x = -x;
        }
        else if ( position == "FFDL" ) // Feet-first decubitus-left (FFDL)
        {
            double v = x;
            x = y;
            y = v;

            y = -y;
        }
        else if ( position == "FFP" ) // Feet-first prone (FFP)
        {
        }
        else if ( position == "FFS" ) // Feet-first supine (FFS)
        {
            x = -x;
            y = -y;
        }
        else 
        {
            GERROR_STREAM("Unknown position string :" << position);
            return false;
        }

        return true;
    }

    bool DeviceCoordinateSystemToPatientCoordinateSystem(double& x, double& y, double& z, const std::string& position)
    {
        if ( position == "HFS" ) // Head-first supine (HFS)
        {
            y = -y;
            z = -z;
        }
        else if ( position == "HFP" ) // Head-first prone (HFP)
        {
            x = -x;
            z = -z;
        }
        else if ( position == "HFDR" ) // Head-first decubitus-right 
        {
            double v = x;
            x = y;
            y = v;
            z = -z;
        }
        else if ( position == "HFDL" ) // Head-first decubitus-left (HFDL)
        {
            double v = x;
            x = y;
            y = v;

            x = -x;
            y = -y;
            z = -z;
        }
        else if ( position == "FFDR" ) // Feet-first decubitus-right (FFDR)
        {
            double v = x;
            x = y;
            y = v;

            y = -y;
        }
        else if ( position == "FFDL" ) // Feet-first decubitus-left (FFDL)
        {
            double v = x;
            x = y;
            y = v;

            x = -x;
        }
        else if ( position == "FFP" ) // Feet-first prone (FFP)
        {
        }
        else if ( position == "FFS" ) // Feet-first supine (FFS)
        {
            x = -x;
            y = -y;
        }
        else 
        {
            GERROR_STREAM("Unknown position string :" << position);
            return false;
        }

        return true;
    }

    std::string getKSpaceFilterNameFromType(Gadgetron::ISMRMRDKSPACEFILTER type)
    {
        std::string name;

        if (type == ISMRMRD_FILTER_GAUSSIAN)
        {
            name = "Gaussian";
        }
        else if (type == ISMRMRD_FILTER_HANNING)
        {
            name = "Hanning";
        }
        else if (type == ISMRMRD_FILTER_TAPERED_HANNING)
        {
            name = "TaperedHanning";
        }
        else
        {
            name = "None";
        }

        return name;
    }
}
