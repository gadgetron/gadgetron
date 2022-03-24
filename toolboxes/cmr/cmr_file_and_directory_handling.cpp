/** \file       cmr_file_and_directory_handling.cpp
    \brief      Functionalities to handle file and directory for CMR.
    \author     Hui Xue
*/

#include "cmr_file_and_directory_handling.h"
#include "io/primitives.h"

#include <boost/filesystem.hpp>

#ifdef _WIN32
#include <windows.h>
#include <Shlwapi.h>
#pragma comment(lib, "shlwapi.lib")
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#endif // _WIN32

namespace Gadgetron
{
    void create_folder_with_all_permissions(const std::string& workingdirectory)
    {
        if (!boost::filesystem::exists(workingdirectory))
        {
            boost::filesystem::path workingPath(workingdirectory);
            try
            {
                boost::filesystem::create_directories(workingPath);
            }
            catch (...)
            {
                GERROR("Error creating the working directory.\n");
                GADGET_THROW("Exceptions happened in create_folder_with_all_permissions(...) ... ");
            }

            // set the permission for the folder
#ifdef _WIN32
            try
            {
                boost::filesystem::permissions(workingPath, boost::filesystem::all_all);
            }
            catch (...)
            {
                GERROR("Error changing the permission of the working directory.\n");
                GADGET_THROW("Exceptions happened in create_folder_with_all_permissions(...) ... ");
            }
#else
            // in case an older version of boost is used in non-win system
            // the system call is used
            int res = chmod(workingPath.string().c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH | S_IXOTH);
            if (res != 0)
            {
                GERROR("Error changing the permission of the working directory.\n");
                GADGET_THROW("Exceptions happened in create_folder_with_all_permissions(...) ... ");
            }
#endif // _WIN32
        }
    }

    void clean_items_older_than(const std::string& workingdirectory, double hours)
    {
        try
        {
            // get current time
            std::time_t curr_time_UTC_;

            std::time(&curr_time_UTC_);
            struct tm* currTm = std::gmtime(&curr_time_UTC_);
            curr_time_UTC_ = std::mktime(currTm);

            // list and clean the content in the workingdirectory
            boost::filesystem::path p(workingdirectory);

            if (boost::filesystem::exists(p))
            {
                if (boost::filesystem::is_directory(p))
                {
                    typedef std::vector<boost::filesystem::path> vec;
                    vec v;
                    v.reserve(100);

                    std::copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), back_inserter(v));
                    std::sort(v.begin(), v.end());

                    GDEBUG_STREAM("A total of " << v.size() << " items are found ... ");

                    // if needed, clean the storage first
                    std::string filename;

                    for (vec::const_iterator it(v.begin()); it != v.end(); ++it)
                    {
                        filename = it->string();

                        // find the file creation/modification time
                        std::time_t lastWriteTime = last_write_time(*it);
                        struct tm* lastWriteTm = std::gmtime(&lastWriteTime);
                        lastWriteTime = std::mktime(lastWriteTm);

                        if (std::abs((double)lastWriteTime - (double)curr_time_UTC_) > hours * 3600.0)
                        {
#ifdef _WIN32
                            boost::filesystem::remove(*it);
#else
                            int res = remove(filename.c_str());
                            if (res != 0)
                            {
                                GERROR_STREAM("clean_items_older_than. error removing " << filename);
                            }
#endif // _WIN32
                        }
                    }

                    // update the file list
                    v.clear();
                    std::copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), back_inserter(v));
                    std::sort(v.begin(), v.end());

                    GDEBUG_STREAM("A total of " << v.size() << " items are found after cleaning ... ");
                }
            }
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in clean_items_older_than(...) ... ");
        }
    }

    void list_items_from_session(const std::string& workingdirectory, const std::string& session_id, std::vector<std::string>& items_list)
    {
        try
        {
            std::time_t curr_time_UTC_;

            std::time(&curr_time_UTC_);
            struct tm* currTm = std::gmtime(&curr_time_UTC_);
            curr_time_UTC_ = std::mktime(currTm);

            boost::filesystem::path p(workingdirectory);

            if (boost::filesystem::exists(p))
            {
                if (boost::filesystem::is_directory(p))
                {
                    typedef std::vector<boost::filesystem::path> vec;
                    vec v;
                    v.reserve(100);

                    std::copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), back_inserter(v));
                    std::sort(v.begin(), v.end());

                    std::vector<std::string> filenames;
                    std::vector<double> changed_time;

                    for (vec::const_iterator it(v.begin()); it != v.end(); ++it)
                    {
                        std::string filename = it->string();

                        if (filename.find(session_id) != std::string::npos)
                        {
                            filenames.push_back(filename);

                            // find the file creation/modification time
                            std::time_t lastWriteTime = last_write_time(*it);
                            struct tm* lastWriteTm = std::gmtime(&lastWriteTime);
                            lastWriteTime = std::mktime(lastWriteTm);

                            changed_time.push_back(std::abs((double)lastWriteTime - (double)curr_time_UTC_));
                        }
                    }

                    // sort the list with new files first
                    std::vector<double> changed_time_sorted(changed_time);
                    std::sort(changed_time_sorted.begin(), changed_time_sorted.end());

                    size_t N = filenames.size();

                    for(size_t ii=0; ii<N; ii++)
                    {
                        size_t jj;
                        for (jj=0; jj<N; jj++)
                        {
                            if(std::abs(changed_time_sorted[ii]- changed_time[jj])<0.1)
                            {
                                break;
                            }
                        }

                        items_list.push_back(filenames[jj]);
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Exceptions happened in list_items_from_session(...) ... ");
        }
    }

    template<typename T>
    void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<T> >& items)
    {
        try
        {
            boost::filesystem::path p(workingdirectory);

            if (boost::filesystem::exists(p))
            {
                if (boost::filesystem::is_directory(p))
                {
                    size_t N = items_list.size();
                    size_t ii;

                    items.clear();

                    if(N==0)
                    {
                        return;
                    }

                    items.resize(N);

                    for (ii = 0; ii < N; ii++)
                    {
                        std::string file_name = workingdirectory + "/" + items_list[ii];

                        if(items_list[ii].find(workingdirectory)!=std::string::npos)
                        {
                            file_name = items_list[ii];
                        }

                        std::ifstream fid;
                        fid.open(file_name.c_str(), std::ios::binary);

                        if (!fid)
                        {
                            GERROR_STREAM("load_items cannot open file stream : " << file_name);
                            continue;
                        }

                        Gadgetron::Core::IO::read(fid, items[ii]);
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Exceptions happened in load_items(...) ... ");
        }
    }

    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDArray<float> > >& items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDArray<double> > >& items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDArray< std::complex<float> > > >& items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDArray< std::complex<double> > > >& items);

    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDImage<float, 2> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDImage<float, 3> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDImage<double, 2> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDImage<double, 3> > > & items);

    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDImage<std::complex<float>, 2> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDImage<std::complex<float>, 3> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDImage<std::complex<double>, 2> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoNDImage<std::complex<double>, 3> > > & items);

    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoMRImage<float, 2> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoMRImage<float, 3> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoMRImage<double, 2> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoMRImage<double, 3> > > & items);

    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoMRImage<std::complex<float>, 2> > >& items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoMRImage<std::complex<float>, 3> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoMRImage<std::complex<double>, 2> > > & items);
    template EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<hoMRImage<std::complex<double>, 3> > > & items);

    template<typename T>
    void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<T>& item)
    {
        try
        {
            boost::filesystem::path p(workingdirectory);

            if (boost::filesystem::exists(p))
            {
                if (boost::filesystem::is_directory(p))
                {
                    std::string file_name = workingdirectory + "/" + session_id;

                    std::ofstream fid;
                    fid.open(file_name.c_str(), std::ios::binary);

                    if (!fid)
                    {
                        GERROR_STREAM("save_item cannot open file stream : " << file_name);
                        GADGET_THROW("Exceptions happened in save_item(...) ... ");
                    }

                    Gadgetron::Core::IO::write(fid, item);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Exceptions happened in save_item(...) ... ");
        }
    }

    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDArray<float> >& item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDArray<double> >& item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDArray< std::complex<float> > >& item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDArray< std::complex<double> > >& item);

    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDImage<float, 2> > & item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDImage<float, 3> > & item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDImage<double, 2> > & item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDImage<double, 3> > & item);

    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDImage<std::complex<float>, 2> > & item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDImage<std::complex<float>, 3> > & item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDImage<std::complex<double>, 2> > & item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoNDImage<std::complex<double>, 3> > & item);

    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoMRImage<float, 2> >& item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoMRImage<float, 3> >& item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoMRImage<double, 2> >& item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoMRImage<double, 3> >& item);

    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoMRImage<std::complex<float>, 2> >& item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoMRImage<std::complex<float>, 3> >& item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoMRImage<std::complex<double>, 2> >& item);
    template EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<hoMRImage<std::complex<double>, 3> >& item);
}
