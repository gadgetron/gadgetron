
/** \file   mri_core_dependencies.cpp
    \brief  Implementation useful utility functionalities for MRI dependency data handling
    \author Hui Xue
*/

#include "mri_core_dependencies.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_utils.h"
#include <fstream>
#include "io/primitives.h"

namespace Gadgetron
{
    template <typename T> 
    void save_dependency_data(const std::string& ismrmrd_header, const hoNDArray<T>& scc_array, const ISMRMRD::AcquisitionHeader& scc_header, 
                            const hoNDArray<T>& body_array, const ISMRMRD::AcquisitionHeader& body_header, const std::string& filename)
    {
        try
        {
            std::ofstream outfile;
            outfile.open( filename.c_str(), std::ios::out | std::ios::binary );

            if (outfile.good())
            {
                GDEBUG_STREAM( "Write out the dependency data file : " << filename );

                size_t len_ismrmrd_header = ismrmrd_header.length();
                outfile.write( reinterpret_cast<char*>(&len_ismrmrd_header), sizeof(size_t) );
                outfile.write( ismrmrd_header.c_str(), len_ismrmrd_header );

                auto array_byte_length = [](auto& array) {return array.get_number_of_bytes()+ sizeof(size_t)*(array.get_number_of_dimensions()+1);};
                size_t len_scc_array = array_byte_length(scc_array);
                size_t len_body_array = array_byte_length(body_array);

                outfile.write( reinterpret_cast<char*>(&len_scc_array), sizeof(size_t) );
                Core::IO::write(outfile,scc_array);

                outfile.write( reinterpret_cast<char*>(&len_body_array), sizeof(size_t) );
                Core::IO::write(outfile,body_array);

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
                scc_array = Core::IO::read<hoNDArray<T>>(infile);

                // -----------------------------

                size_t len_body_array = 0;
                infile.read( reinterpret_cast<char*>(&len_body_array), sizeof(size_t));
                body_array = Core::IO::read<hoNDArray<T>>(infile);


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
            GADGET_THROW("Errors in load_dependency_data(...) ... ");
        }
    }

    template EXPORTMRICORE void load_dependency_data(const std::string& filename, std::string& ismrmd_header, hoNDArray< std::complex<float> >& scc_array, ISMRMRD::AcquisitionHeader& scc_header, hoNDArray< std::complex<float> >& body_array, ISMRMRD::AcquisitionHeader& body_header);
    template EXPORTMRICORE void load_dependency_data(const std::string& filename, std::string& ismrmd_header, hoNDArray< std::complex<double> >& scc_array, ISMRMRD::AcquisitionHeader& scc_header, hoNDArray< std::complex<double> >& body_array, ISMRMRD::AcquisitionHeader& body_header);
}
