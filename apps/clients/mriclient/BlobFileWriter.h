#ifndef BLOB_FILE_WRITER_H
#define BLOB_FILE_WRITER_H

#include <fstream>
#include <iomanip>

#include "GadgetMessageInterface.h"

namespace Gadgetron {

#define MAX_BLOBS_LOG_10    6

class BlobFileWriter : public GadgetMessageReader
{

    public:
        BlobFileWriter(std::string fileprefix, std::string filesuffix)
            : number_of_calls_(0)
            , file_prefix(fileprefix)
            , file_suffix(filesuffix)
        {
        }

        virtual ~BlobFileWriter() {};

        virtual ACE_Message_Block* read(ACE_SOCK_Stream* socket)
        {
            ssize_t recv_count = 0;

            // MUST READ 32-bits
            uint32_t nbytes;
            if ((recv_count = socket->recv_n(&nbytes, sizeof(nbytes))) <= 0) {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, BlobFileWriter, failed to read Blob Header\n")) );
                return 0;
            }

            char *data = new char[nbytes];
            if ((recv_count = socket->recv_n(data, nbytes)) <= 0) {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, BlobFileWriter, failed to read blob from socket\n")) );
                return 0;
            }

            if (this->process_image(nbytes, data) < 0) {
                GADGET_DEBUG1("Failed to process image\n");
                return 0;
            }

            delete[] data;

            // The GadgetronConnector expects an ACE_Message_Block* (NOT NULL)
            ACE_Message_Block *mb = new ACE_Message_Block();

            return mb;
        }

        virtual int process_image(const unsigned int bytes, const char* data)
        {
            std::stringstream filename;

            // Create the filename: (prefix_%06.suffix)
            filename << file_prefix << "_";
            filename << std::setfill('0') << std::setw(MAX_BLOBS_LOG_10) << number_of_calls_;
            filename << "." << file_suffix;

            std::ofstream outfile;
            outfile.open (filename.str().c_str(), std::ios::out|std::ios::binary);

            ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Writing image %s\n"), filename.str().c_str()) );

            if (outfile.good()) {
                /* write 'size' bytes starting at 'data's pointer */
                outfile.write(data, bytes);
                outfile.close();
                number_of_calls_++;
            } else {
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

#endif //BLOB_FILE_WRITER_H
