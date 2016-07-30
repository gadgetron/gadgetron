
/** \file   mri_core_dependencies.cpp
    \brief  Implementation useful utility functionalities for MRI dependency data handling
    \author Hui Xue
*/

#include "mri_core_dependencies.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_utils.h"
#include <fstream>

namespace Gadgetron
{
    template <typename T> 
    void save_dependency_data(const std::string& ismrmrd_header, const hoNDArray<T>& scc_array, const ISMRMRD::AcquisitionHeader& scc_header, 
                            const hoNDArray<T>& body_array, const ISMRMRD::AcquisitionHeader& body_header, const std::string& filename)
    {
        char* buf = NULL;
        char* buf_scc_array = NULL;
        char* buf_body_array = NULL;

        try
        {
            std::ofstream outfile;
            outfile.open( filename.c_str(), std::ios::out | std::ios::binary );

            if (outfile.good())
            {
                GDEBUG_STREAM( "Write out the dependency data file : " << filename );

                size_t len_ismrmrd_header = ismrmrd_header.length();

                size_t len_scc_array = 0;
                GADGET_CHECK_THROW(scc_array.serialize(buf_scc_array, len_scc_array));

                size_t len_body_array = 0;
                GADGET_CHECK_THROW(scc_array.serialize(buf_body_array, len_body_array));

                outfile.write( reinterpret_cast<char*>(&len_ismrmrd_header), sizeof(size_t) );
                outfile.write( ismrmrd_header.c_str(), len_ismrmrd_header );

                outfile.write( reinterpret_cast<char*>(&len_scc_array), sizeof(size_t) );
                outfile.write( buf_scc_array, len_scc_array );

                outfile.write( reinterpret_cast<char*>(&len_body_array), sizeof(size_t) );
                outfile.write( buf_body_array, len_body_array );

                outfile.write( reinterpret_cast<const char*>(&scc_header), sizeof(ISMRMRD::AcquisitionHeader) );
                outfile.write( reinterpret_cast<const char*>(&body_header), sizeof(ISMRMRD::AcquisitionHeader) );

                outfile.close();
            }
            else
            {
                GERROR_STREAM("Failed to open dependency data file for writing : " << filename);
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in save_dependency_data(...) ... ");

            if(buf!=NULL) delete [] buf;
            if(buf_scc_array!=NULL) delete [] buf_scc_array;
            if(buf_body_array!=NULL) delete [] buf_body_array;
        }
    }

    template EXPORTMRICORE void save_dependency_data(const std::string& ismrmd_header, const hoNDArray< std::complex<float> >& scc_array, const ISMRMRD::AcquisitionHeader& scc_header, const hoNDArray< std::complex<float> >& body_array, const ISMRMRD::AcquisitionHeader& body_header, const std::string& filename);
    template EXPORTMRICORE void save_dependency_data(const std::string& ismrmd_header, const hoNDArray< std::complex<double> >& scc_array, const ISMRMRD::AcquisitionHeader& scc_header, const hoNDArray< std::complex<double> >& body_array, const ISMRMRD::AcquisitionHeader& body_header, const std::string& filename);

    // ------------------------------------------------------------------------

    template <typename T> 
    void load_dependency_data(const std::string& filename, std::string& ismrmd_header, hoNDArray<T>& scc_array, ISMRMRD::AcquisitionHeader& scc_header, 
                        hoNDArray<T>& body_array, ISMRMRD::AcquisitionHeader& body_header)
    {
        try
        {
            std::ifstream infile;
            infile.open (filename.c_str(), std::ios::in|std::ios::binary);

            if (infile.good())
            {
                size_t xml_length;
                infile.read( reinterpret_cast<char*>(&xml_length), sizeof(size_t));

                std::string xml_str(xml_length,'\0');
                infile.read(const_cast<char*>(xml_str.c_str()), xml_length);
                ismrmd_header = xml_str;

                // -----------------------------

                size_t len_scc_array = 0;
                infile.read( reinterpret_cast<char*>(&len_scc_array), sizeof(size_t));

                std::vector<char> buf_scc(len_scc_array);
                infile.read( &buf_scc[0], len_scc_array);

                scc_array.deserialize(&buf_scc[0], len_scc_array);

                // -----------------------------

                size_t len_body_array = 0;
                infile.read( reinterpret_cast<char*>(&len_body_array), sizeof(size_t));

                std::vector<char> buf_body(len_body_array);
                infile.read( &buf_body[0], len_body_array);

                body_array.deserialize(&buf_scc[0], len_scc_array);

                // -----------------------------

                infile.read( reinterpret_cast<char*>(&scc_header), sizeof(ISMRMRD::AcquisitionHeader));
                infile.read( reinterpret_cast<char*>(&body_header), sizeof(ISMRMRD::AcquisitionHeader));

                infile.close();
            }
            else
            {
                GERROR_STREAM("Failed to open dependency data file for reading : " << filename);
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in save_dependency_data(...) ... ");
        }
    }

    template EXPORTMRICORE void load_dependency_data(const std::string& filename, std::string& ismrmd_header, hoNDArray< std::complex<float> >& scc_array, ISMRMRD::AcquisitionHeader& scc_header, hoNDArray< std::complex<float> >& body_array, ISMRMRD::AcquisitionHeader& body_header);
    template EXPORTMRICORE void load_dependency_data(const std::string& filename, std::string& ismrmd_header, hoNDArray< std::complex<double> >& scc_array, ISMRMRD::AcquisitionHeader& scc_header, hoNDArray< std::complex<double> >& body_array, ISMRMRD::AcquisitionHeader& body_header);
}
