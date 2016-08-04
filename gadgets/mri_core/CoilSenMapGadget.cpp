
#include "CoilSenMapGadget.h"
#include "hoArmadillo.h"
#include "hoNDArray_elemwise.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_reductions.h"

#include "mri_core_dependencies.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

#ifndef _WIN32
    #include <sys/types.h>
    #include <sys/stat.h>
#endif // _WIN32

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

namespace Gadgetron {

    CoilSenMapGadget::CoilSenMapGadget()
    {
    }

    CoilSenMapGadget::~CoilSenMapGadget()
    {

    }

    int CoilSenMapGadget::process_config( ACE_Message_Block* mb )
    {
        if (!debug_folder.value().empty())
        {
            Gadgetron::get_debug_folder_path(debug_folder.value(), debug_folder_full_path_);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is " << debug_folder_full_path_);
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is not set ... ");
        }

        if (!workingDirectory.value().empty())
        {
            coil_sen_dependency_folder_ = workingDirectory.value();
        }
        else
        {
#ifdef _WIN32
            coil_sen_dependency_folder_ = std::string( "c:\\temp\\gadgetron\\" );
#else
            coil_sen_dependency_folder_ = std::string( "/tmp/gadgetron/" );
#endif // _WIN32
        }

        GDEBUG( "Folder to store coil map dependencies is %s\n", coil_sen_dependency_folder_.c_str() );

        // deserialize the ismrmrd header
        ISMRMRD::deserialize(mb->rd_ptr(), current_ismrmrd_header_);

        if (current_ismrmrd_header_.measurementInformation)
        {
            if (current_ismrmrd_header_.measurementInformation->measurementID)
            {
                measurement_id_ = *current_ismrmrd_header_.measurementInformation->measurementID;
                GDEBUG( "Measurement ID is %s\n", measurement_id_.c_str() );
            }
            else
            {
                GERROR_STREAM("Incoming ismrmrd header does not have measurementInformation.measurementID ... ");
                return GADGET_FAIL;
            }
        }
        else
        {
            GERROR_STREAM("Incoming ismrmrd header does not have measurementInformation ... ");
            return GADGET_FAIL;
        }

        full_name_stored_coil_sen_dependency_ = this->generate_coil_sen_dependency_filename( measurement_id_ );
        GDEBUG( "Stored coil sen dependency is %s\n", full_name_stored_coil_sen_dependency_.c_str() );

        //// find the measurementID of this scan
        //if (current_ismrmrd_header_.measurementInformation)
        //{
        //    if (current_ismrmrd_header_.measurementInformation->measurementID)
        //    {
        //        measurement_id_ = *current_ismrmrd_header_.measurementInformation->measurementID;
        //        GDEBUG( "Measurement ID is %s\n", measurement_id_.c_str() );
        //    }

        //    // find the noise depencies if any
        //    if (current_ismrmrd_header_.measurementInformation->measurementDependency.size() > 0)
        //    {
        //        measurement_id_of_coil_sen_dependency_.clear();

        //        std::vector<ISMRMRD::MeasurementDependency>::const_iterator iter = current_ismrmrd_header_.measurementInformation->measurementDependency.begin();
        //        for (; iter != current_ismrmrd_header_.measurementInformation->measurementDependency.end(); iter++)
        //        {
        //            std::string dependencyType = iter->dependencyType;
        //            std::string dependencyID = iter->measurementID;

        //            GDEBUG( "Found dependency measurement : %s with ID %s\n", dependencyType.c_str(), dependencyID.c_str() );

        //            if (dependencyType == "SenMap" || dependencyType == "senmap") {
        //                measurement_id_of_coil_sen_dependency_ = dependencyID;
        //            }
        //        }

        //        if (!measurement_id_of_coil_sen_dependency_.empty())
        //        {
        //            GDEBUG( "Measurement ID of coil sen dependency is %s\n", measurement_id_of_coil_sen_dependency_.c_str() );

        //            full_name_stored_coil_sen_dependency_ = this->generate_coil_sen_dependency_filename( measurement_id_of_coil_sen_dependency_ );
        //            GDEBUG( "Stored coil sen dependency is %s\n", full_name_stored_coil_sen_dependency_.c_str() );
        //        }
        //    }
        //}

        return GADGET_OK;
    }

    std::string CoilSenMapGadget::generate_coil_sen_dependency_filename( const std::string& measurement_id )
    {
        std::string full_name_stored_coil_sen_dependency;

        full_name_stored_coil_sen_dependency = coil_sen_dependency_folder_;
        full_name_stored_coil_sen_dependency.append( "/" );
        full_name_stored_coil_sen_dependency.append( coil_sen_dependency_prefix.value() );
        full_name_stored_coil_sen_dependency.append( "_" );
        full_name_stored_coil_sen_dependency.append( measurement_id );

        return full_name_stored_coil_sen_dependency;
    }

    bool CoilSenMapGadget::save_coil_sen_dependency()
    {
        // ismrmrd header
        std::stringstream xml_ss;
        ISMRMRD::serialize( current_ismrmrd_header_, xml_ss );
        std::string xml_str = xml_ss.str();
        uint32_t xml_length = static_cast<uint32_t>(xml_str.size());

        if(body_header_.empty())
        {
            ISMRMRD::AcquisitionHeader header;
            memset(&header, 0, sizeof(ISMRMRD::AcquisitionHeader));

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::save_dependency_data(xml_str, scc_array_, scc_header_[0], body_array_, header, full_name_stored_coil_sen_dependency_));
        }
        else
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::save_dependency_data(xml_str, scc_array_, scc_header_[0], body_array_, body_header_[0], full_name_stored_coil_sen_dependency_));
        }

        //std::string xml_loaded;
        //ISMRMRD::AcquisitionHeader hs, hb;
        //hoNDArray< std::complex<float> > sa, ba;
        //Gadgetron::load_dependency_data(full_name_stored_coil_sen_dependency_, xml_loaded, sa, hs, ba, hb);

        //if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(sa, debug_folder_full_path_ + "CoilSenMap_sa_array");
        //if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(ba, debug_folder_full_path_ + "CoilSenMap_ba_array");

        // set the permission
#ifndef _WIN32
        int res = chmod( full_name_stored_coil_sen_dependency_.c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH | S_IXOTH );
        if (res != 0)
        {
            GDEBUG( "Changing coil sen file permission failed ...\n" );
        }
#endif // _WIN32

        return true;
    }

    int CoilSenMapGadget::process( GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1 )
    {
        // find the dimension limit of stored data

        size_t RO(0), E1(0), E2(0), SLC(0), PHS(0), CON(0), REP(0), SET(0), SEG(0), AVE(0);

        IsmrmrdAcquisitionBucket& bucket = *m1->getObjectPtr();

        size_t N = bucket.data_.size();

        RO = bucket.data_[0].head_->getObjectPtr()->number_of_samples;
        E1 = *std::max_element(bucket.datastats_[0].kspace_encode_step_1.begin(), bucket.datastats_[0].kspace_encode_step_1.end()) + 1;
        E2 = *std::max_element(bucket.datastats_[0].kspace_encode_step_2.begin(), bucket.datastats_[0].kspace_encode_step_2.end()) + 1;
        SLC = *std::max_element(bucket.datastats_[0].slice.begin(), bucket.datastats_[0].slice.end()) + 1;
        PHS = *std::max_element(bucket.datastats_[0].phase.begin(), bucket.datastats_[0].phase.end()) + 1;
        CON = *std::max_element(bucket.datastats_[0].contrast.begin(), bucket.datastats_[0].contrast.end()) + 1;
        REP = *std::max_element(bucket.datastats_[0].repetition.begin(), bucket.datastats_[0].repetition.end()) + 1;
        SET = *std::max_element(bucket.datastats_[0].set.begin(), bucket.datastats_[0].set.end()) + 1;
        SEG = *std::max_element(bucket.datastats_[0].segment.begin(), bucket.datastats_[0].segment.end()) + 1;
        AVE = *std::max_element(bucket.datastats_[0].average.begin(), bucket.datastats_[0].average.end()) + 1;

        GDEBUG_CONDITION_STREAM(verbose.value(), "CoilSenMapGadget, received buckert [RO E1 E2 SLC PHS CON REP SET SEG AVE] : [" 
                                                                        << RO << " " << E1 << " " << E2 << " " 
                                                                        << SLC << " " << PHS << " " << CON << " " 
                                                                        << REP << " " << SET << " " << SEG << " " << SLC << "]");

        // find max CHA for each SET
        size_t CHA_1st_set(0), CHA_2nd_set(0);

        size_t n;
        for (n=0; n<N; n++)
        {
            if(bucket.data_[n].head_->getObjectPtr()->idx.set==0)
            {
                if(CHA_1st_set < bucket.data_[n].head_->getObjectPtr()->active_channels)
                {
                    CHA_1st_set = bucket.data_[n].head_->getObjectPtr()->active_channels;
                }
            }

            if(bucket.data_[n].head_->getObjectPtr()->idx.set==1)
            {
                if(CHA_2nd_set < bucket.data_[n].head_->getObjectPtr()->active_channels)
                {
                    CHA_2nd_set = bucket.data_[n].head_->getObjectPtr()->active_channels;
                }
            }
        }

        GDEBUG_CONDITION_STREAM(verbose.value(), "CoilSenMapGadget, received number of channel, set 0 and set 1 : " << CHA_1st_set << " and " << CHA_2nd_set);

        // allocate surface coil and body coil array
        std::vector<size_t> dim(11, 1);
        dim[0] = RO;
        dim[1] = E1;
        dim[2] = E2;
        dim[3] = CHA_1st_set;
        dim[4] = SLC;
        dim[5] = PHS;
        dim[6] = CON;
        dim[7] = REP;
        dim[8] = 1;
        dim[9] = 1;
        dim[10] = AVE;

        scc_array_.create(dim);
        Gadgetron::clear(scc_array_);

        if(CHA_2nd_set>0)
        {
            std::vector<size_t> dim_body(dim);
            dim_body[3] = CHA_2nd_set;

            body_array_.create(dim_body);
            Gadgetron::clear(body_array_);
        }

        // fill the array

        std::vector<size_t> ind(11);

        size_t cha;
        for (n=0; n<N; n++)
        {
            ind[0] = 0;
            ind[1] = bucket.data_[n].head_->getObjectPtr()->idx.kspace_encode_step_1;
            ind[2] = bucket.data_[n].head_->getObjectPtr()->idx.kspace_encode_step_2;

            ind[4] = bucket.data_[n].head_->getObjectPtr()->idx.slice;
            ind[5] = bucket.data_[n].head_->getObjectPtr()->idx.phase;
            ind[6] = bucket.data_[n].head_->getObjectPtr()->idx.contrast;
            ind[7] = bucket.data_[n].head_->getObjectPtr()->idx.repetition;
            ind[8] = 0;
            ind[9] = 0;
            ind[10] = bucket.data_[n].head_->getObjectPtr()->idx.average;

            if(bucket.data_[n].head_->getObjectPtr()->idx.set==0)
            {
                scc_header_.push_back(*bucket.data_[n].head_->getObjectPtr());

                for (cha=0; cha<CHA_1st_set; cha++)
                {
                    ind[3] = cha;
                    size_t array_index = scc_array_.calculate_offset(ind);
                    std::complex<float>* pData = &(scc_array_(array_index));

                    memcpy(pData, bucket.data_[n].data_->getObjectPtr()->begin() + cha*RO, sizeof(std::complex<float>)*RO);
                }
            }

            if(bucket.data_[n].head_->getObjectPtr()->idx.set==1)
            {
                body_header_.push_back(*bucket.data_[n].head_->getObjectPtr());

                for (cha=0; cha<CHA_2nd_set; cha++)
                {
                    ind[3] = cha;
                    size_t array_index = body_array_.calculate_offset(ind);
                    std::complex<float>* pData = &(body_array_(array_index));

                    memcpy(pData, bucket.data_[n].data_->getObjectPtr()->begin() + cha*RO, sizeof(std::complex<float>)*RO);
                }
            }
        }

        if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(scc_array_, debug_folder_full_path_ + "CoilSenMap_scc_array");

        if(CHA_2nd_set>0)
        {
            if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(body_array_, debug_folder_full_path_ + "CoilSenMap_body_array");
        }

        if (!this->save_coil_sen_dependency())
        {
            GERROR_STREAM("CoilSenMapGadget, save_coil_sen_dependency failed ... ");
        }

        m1->release();

        return GADGET_OK;

    }

    int CoilSenMapGadget::close( unsigned long flags )
    {
        GDEBUG_CONDITION_STREAM(true, "CoilSenMapGadget - close(flags) : " << flags);

        if (BaseClass::close( flags ) != GADGET_OK) return GADGET_FAIL;

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE( CoilSenMapGadget )

} // namespace Gadgetron
