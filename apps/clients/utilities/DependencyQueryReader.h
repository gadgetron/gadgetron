
/** \file   DependencyQueryReader.h
    \brief  Implement the writer to write the dependency query reults into a file
    \author Hui Xue
*/

#pragma once

#include <fstream>
#include <iomanip>

#include "GadgetMessageInterface.h"
#include "ismrmrd/meta.h"

namespace Gadgetron
{

class DependencyQueryReader : public GadgetMessageReader
{
    public:

        DependencyQueryReader(std::string filename) : number_of_calls_(0) , filename_(filename)
        {
        }

        virtual ~DependencyQueryReader()
        {
        }

        virtual ACE_Message_Block* read(ACE_SOCK_Stream* socket)
        {
            ssize_t recv_count = 0;

            typedef unsigned long long size_t_type;

            size_t_type len(0);
            if ( ( recv_count = socket->recv_n( &len, sizeof(size_t_type)) ) <= 0 )
            {
	      GERROR("DependencyQueryReader, failed to read query results length\n");
	      return 0;
            }

            char* buf = NULL;
            try
            {
                buf = new char[len];
                if ( buf == NULL )
                {
		  GERROR("DependencyQueryReader, failed to allocate buffer\n");
		  return 0;
                }

                memset(buf, '\0', len);
                memcpy(buf, &len, sizeof(size_t_type));
            }
            catch (std::runtime_error &err)
            {
                GEXCEPTION(err,"DependencyQueryReader, failed to allocate buffer\n");
                return 0;
            }

            if ( ( recv_count = socket->recv_n( buf, len) ) <= 0 )
            {
	      GERROR("DependencyQueryReader, failed to read query results\n");
	      delete [] buf;
	      return 0;
            }

            std::ofstream outfile;
            outfile.open (filename_.c_str(), std::ios::out|std::ios::binary);

            if (outfile.good())
            {
                outfile.write(buf, len);
                outfile.close();
                number_of_calls_++;
            }
            else
            {
                delete[] buf;

                GERROR_STREAM("File " << filename_ << " is not good for writing\n");
                return 0;
            }

            delete[] buf;

            // The GadgetronConnector expects an ACE_Message_Block* (NOT NULL)
            ACE_Message_Block *mb = new ACE_Message_Block();

            return mb;
        }

    protected:

        size_t number_of_calls_;
        std::string filename_;
};

} // namespace Gadgetron
