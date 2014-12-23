#include "Gadgetron.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "DependencyQueryGadget.h"
#include "gtPlusUtil.h"

#include <boost/version.hpp>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;

namespace Gadgetron
{
    DependencyQueryGadget::DependencyQueryGadget()
    {
        processed_in_close_ = false;

        noise_dependency_prefix_ = "GadgetronNoiseCovarianceMatrix";

        noise_dependency_attrib_name_ = "NoiseDependencies";

        clean_storage_while_query_ = true;
        time_limit_in_storage_ = 24.0;

        // get current time
        std::time(&curr_time_UTC_);
        struct tm* currTm = std::gmtime(&curr_time_UTC_);
        curr_time_UTC_ = std::mktime(currTm);
    }

    DependencyQueryGadget::~DependencyQueryGadget()
    {
    }

    int DependencyQueryGadget::process_config(ACE_Message_Block* mb)
    {
        return GADGET_OK;
    }

    int DependencyQueryGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< ValueType > >* m2)
    {
        return GADGET_OK;
    }

    int DependencyQueryGadget::close(unsigned long flags)
    {
        typedef unsigned long long size_t_type;

        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

        if ( !processed_in_close_ )
        {
            processed_in_close_ = true;

            boost::shared_ptr<std::string> str = this->get_string_value("workingDirectory");
            if ( !str->empty() )
            {
                noise_dependency_folder_ = *str;
            }
            else
            {
	      //This is an error, we should not be writing dependencies without having a working directory
	      return GADGET_FAIL;
            }
            GADGET_MSG("Folder to store noise dependencies is " << noise_dependency_folder_);

            str = this->get_string_value("noise_dependency_prefix");

            if ( !str->empty() )
            {
                noise_dependency_prefix_ = *str;
            }

            str = this->get_string_value("noise_dependency_attrib_name");

            if ( !str->empty() )
            {
                noise_dependency_attrib_name_ = *str;
            }

            clean_storage_while_query_ = this->get_bool_value("clean_storage_while_query");
            GADGET_MSG( "clean_storage_while_query_ is " << clean_storage_while_query_);

            time_limit_in_storage_ = this->get_double_value("time_limit_in_storage");
            if ( time_limit_in_storage_ < 0 )
            {
                time_limit_in_storage_ = 24.0;
            }
            GADGET_MSG( "time_limit_in_storage_ is " << time_limit_in_storage_);

            // list the content in the noise dependency folder
            path p (noise_dependency_folder_);

            try
            {
                if ( exists(p) )
                {
                    if ( is_directory(p) )
                    {
                        typedef std::vector<path> vec;
                        vec v;
                        v.reserve(100);

                        copy(directory_iterator(p), directory_iterator(), back_inserter(v));
                        sort(v.begin(), v.end());

                        GADGET_MSG( "A total of " << v.size() << " dependency measurements are found ... ");

                        // if needed, clean the storage first
                        std::string filename;

                        if ( clean_storage_while_query_ )
                        {
                            Gadgetron::gtPlus::gtPlusUtil<ValueType> gt_util;

                            for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
                            {
                                filename = it->string();

                                // find the file creation/modification time
                                std::time_t lastWriteTime = last_write_time(*it);
                                struct tm* lastWriteTm = std::gmtime(&lastWriteTime);
                                lastWriteTime = std::mktime(lastWriteTm);

                                if ( GT_ABS( (double)lastWriteTime - (double)curr_time_UTC_ ) > time_limit_in_storage_*3600.0 )
                                {
                                    remove(*it);
                                }
                            }

                            // update the file list
                            v.clear();
                            copy(directory_iterator(p), directory_iterator(), back_inserter(v));
                            sort(v.begin(), v.end());

                            GADGET_MSG( "A total of " << v.size() << " dependency measurements are found after cleaning ... ");
                        }

                        // declear the attributes
                        Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>* m1 = new Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>();

                        size_t count = 0;
                        size_t ind;

                        for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
                        {
#                       if BOOST_VERSION < 104600
                            filename = it->filename();
#                       else
                            filename = it->filename().string();
#                       endif
                            ind = filename.find(noise_dependency_prefix_);

                            if ( ind != std::string::npos )
                            {
                                m1->getObjectPtr()->append(noise_dependency_attrib_name_.c_str(), filename.c_str());
                                count++;
                            }
                        }

                        GADGET_MSG( "A total of " << count << " noise dependency measurements are found ... ");

                        if ( count == 0 )
                        {
                            // put into a dummy item
                            m1->getObjectPtr()->set(noise_dependency_attrib_name_.c_str(), "Dummy");
                        }

                        // send the found dependencies
                        GadgetContainerMessage<GadgetMessageIdentifier>* mb = new GadgetContainerMessage<GadgetMessageIdentifier>();
                        mb->getObjectPtr()->id = GADGET_MESSAGE_DEPENDENCY_QUERY;
                        mb->cont(m1);

                        int ret =  this->controller_->output_ready(mb);
                        if ( (ret < 0) )
                        {
                            GADGET_DEBUG1("Failed to return massage to controller\n");
                            return GADGET_FAIL;
                        }
                    }
                    else
                    {
                        GADGET_ERROR_MSG( noise_dependency_folder_ << " is not a valid folder ... ");
                    }
                }
                else
                {
                    GADGET_ERROR_MSG("Cannot find dependency folder : " << noise_dependency_folder_);
                }
            }
            catch (const filesystem_error& ex)
            {
                GADGET_ERROR_MSG( ex.what() );
            }
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(DependencyQueryGadget)

} // namespace Gadgetron
