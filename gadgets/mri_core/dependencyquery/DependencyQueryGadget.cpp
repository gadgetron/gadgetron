#include "DependencyQueryGadget.h"

#include <boost/version.hpp>
#include <boost/filesystem.hpp>
#include "Dependency.h"

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

            if ( !workingDirectory.value().empty() )
            {
	      noise_dependency_folder_ = workingDirectory.value();
            }
            else
            {
	      //This is an error, we should not be writing dependencies without having a working directory
	      return GADGET_FAIL;
            }
            GDEBUG_STREAM("Folder to store noise dependencies is " << noise_dependency_folder_);

            if ( !noise_dependency_prefix.value().empty() )
            {
	      noise_dependency_prefix_ = noise_dependency_prefix.value();
            }

            if ( !noise_dependency_attrib_name.value().empty() )
            {
	      noise_dependency_attrib_name_ = noise_dependency_attrib_name.value();
            }

            clean_storage_while_query_ = clean_storage_while_query.value();
            GDEBUG_STREAM( "clean_storage_while_query_ is " << clean_storage_while_query_);

            time_limit_in_storage_ = time_limit_in_storage.value();
            if ( time_limit_in_storage_ < 0 )
            {
                time_limit_in_storage_ = 24.0;
            }
            GDEBUG_STREAM( "time_limit_in_storage_ is " << time_limit_in_storage_);

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

                        GDEBUG_STREAM( "A total of " << v.size() << " dependency measurements are found ... ");

                        // if needed, clean the storage first

                        if ( clean_storage_while_query_ )
                        {
                            for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
                            {
                                std::string filename = it->string();

                                // find the file creation/modification time
                                std::time_t lastWriteTime = last_write_time(*it);
                                struct tm* lastWriteTm = std::gmtime(&lastWriteTime);
                                lastWriteTime = std::mktime(lastWriteTm);

                                if ( std::abs( (double)lastWriteTime - (double)curr_time_UTC_ ) > time_limit_in_storage_*3600.0 )
                                {
                                    remove(*it);
                                }
                            }

                            // update the file list
                            v.clear();
                            copy(directory_iterator(p), directory_iterator(), back_inserter(v));
                            sort(v.begin(), v.end());

                            GDEBUG_STREAM( "A total of " << v.size() << " dependency measurements are found after cleaning ... ");
                        }

                        // declear the attributes
                        auto message = new GadgetContainerMessage<DependencyQuery::Dependency>();

                        auto& dependencies = message->getObjectPtr()->dependencies;

                        size_t count = 0;
                        size_t ind;

                        for (vec::const_iterator it (v.begin()); it != v.end(); ++it)
                        {
                            std::string filename = it->filename().string();
                            ind = filename.find(noise_dependency_prefix_);

                            if ( ind != std::string::npos )
                            {
                                dependencies.append(noise_dependency_attrib_name_.c_str(), filename.c_str());
                                count++;
                            }
                        }

                        GDEBUG_STREAM( "A total of " << count << " noise dependency measurements are found ... ");

                        if ( count == 0 )
                        {
                            // put into a dummy item
                            dependencies.append(noise_dependency_attrib_name_.c_str(), "Dummy");
                        }

                        this->next()->putq(message);
                    }
                    else
                    {
                        GERROR_STREAM( noise_dependency_folder_ << " is not a valid folder ... ");
                    }
                }
                else
                {
                    GERROR_STREAM("Cannot find dependency folder : " << noise_dependency_folder_);
                }
            }
            catch (const filesystem_error& ex)
            {
                GERROR_STREAM( ex.what() );
            }
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(DependencyQueryGadget)

} // namespace Gadgetron
