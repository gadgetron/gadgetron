
/** \file   BlobFileWithAttribWriter.h
    \brief  Implement the blob file writer with filename and meta attributes
    \author Hui Xue
*/

#pragma once

#include <fstream>
#include <iomanip>

#include "GadgetMessageInterface.h"
#include "ismrmrd_meta.h"

namespace Gadgetron
{

#define MAX_BLOBS_LOG_10    6

class BlobFileWithAttribWriter : public GadgetMessageReader
{

    public:
        BlobFileWithAttribWriter(std::string fileprefix, std::string filesuffix)
            : number_of_calls_(0)
            , file_prefix(fileprefix)
            , file_suffix(filesuffix)
        {
        }

        virtual ~BlobFileWithAttribWriter() {};

        virtual ACE_Message_Block* read(ACE_SOCK_Stream* socket)
        {
            ssize_t recv_count = 0;

            // MUST READ 32-bits
            uint32_t nbytes;
            if ((recv_count = socket->recv_n(&nbytes, sizeof(nbytes))) <= 0)
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, BlobFileWithNameWriter, failed to read Blob Header\n")) );
                return 0;
            }

            char *data = new char[nbytes];
            if ((recv_count = socket->recv_n(data, nbytes)) <= 0)
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, BlobFileWithNameWriter, failed to read blob from socket\n")) );
                return 0;
            }

            unsigned long long fileNameLen;
            if ((recv_count = socket->recv_n(&fileNameLen, sizeof(unsigned long long))) <= 0)
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, BlobFileWithNameWriter, failed to read Blob filename length\n")) );
                return 0;
            }

            char* fileNameBuf = new char[fileNameLen+1];
            memset(fileNameBuf, '\0', fileNameLen+1);
            if ((recv_count = socket->recv_n(fileNameBuf, fileNameLen)) <= 0)
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, BlobFileWithNameWriter, failed to read blob filename frome socket\n")) );
                return 0;
            }

            // meta info
            unsigned long long metaLen;
            if ((recv_count = socket->recv_n(&metaLen, sizeof(unsigned long long))) <= 0)
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, BlobFileWithNameWriter, failed to read Blob meta attribute length\n")) );
                return 0;
            }

            char* metaBuf = new char[metaLen+1];
            memset(metaBuf, '\0', metaLen+1);
            if ((recv_count = socket->recv_n(metaBuf, metaLen-sizeof(unsigned long long))) <= 0)
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, BlobFileWithNameWriter, failed to read blob meta attributes frome socket\n")) );
                return 0;
            }

            std::string fileName(fileNameBuf);
            std::string metaAttrib(metaBuf);
            if (this->process_image(nbytes, data, fileName, metaAttrib) < 0)
            {
                GADGET_DEBUG1("Failed to process image\n");
                return 0;
            }

            delete[] data;
            delete[] fileNameBuf;
            delete[] metaBuf;

            // The GadgetronConnector expects an ACE_Message_Block* (NOT NULL)
            ACE_Message_Block *mb = new ACE_Message_Block();

            return mb;
        }

        virtual int process_image(const unsigned int bytes, const char* data, std::string& fileName, std::string& metaAttrib)
        {
            std::string filename_image, filename_attrib;

            // Create the filename: (prefix_%06.suffix)
            if ( file_prefix.empty() )
            {
                filename_image =  fileName + "." + file_suffix;
                filename_attrib =  fileName + "_attrib.xml";
            }
            else
            {
                filename_image = file_prefix + "_" + fileName + "." + file_suffix;
                filename_attrib = file_prefix + "_" + fileName + "_attrib.xml";
            }

            ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Writing image %s\n"), filename_image.c_str()) );

            std::ofstream outfile;
            outfile.open (filename_image.c_str(), std::ios::out|std::ios::binary);

            std::ofstream outfile_attrib;
            outfile_attrib.open (filename_attrib.c_str(), std::ios::out|std::ios::binary);

            if (outfile.good())
            {
                /* write 'size' bytes starting at 'data's pointer */
                outfile.write(data, bytes);
                outfile.close();

                outfile_attrib.write(metaAttrib.c_str(), metaAttrib.length());
                outfile_attrib.close();

                number_of_calls_++;
            }
            else
            {
                GADGET_DEBUG1("File is not good for writing\n");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

    protected:
        size_t number_of_calls_;
        std::string file_prefix;
        std::string file_suffix;
};

} // namespace Gadgetron
