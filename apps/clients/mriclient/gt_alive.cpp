#include "ace/Log_Msg.h"
#include "ace/Get_Opt.h"
#include "ace/OS_NS_string.h"

#include "GadgetronConnector.h"
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "ImageWriter.h"
#include "HDF5ImageWriter.h"
#include "FileInfo.h"
#include "ismrmrd_hdf5.h"
#include "GadgetIsmrmrdReadWrite.h"

#include <fstream>
#include <time.h>
#include <iomanip>

using namespace Gadgetron;


int ACE_TMAIN(int argc, ACE_TCHAR *argv[] )
{
	GadgetronConnector con;

	std::string host("localhost");
	std::string port("9002");

	if (argc > 1) {
		host = std::string(argv[1]);
	}

	if (argc > 2) {
		port = std::string(argv[2]);
	}

	if (con.open(host,port) != 0) {
		ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to connect to the Gadgetron host")));
		return -1;
	}

	//Tell Gadgetron which XML configuration to run.
	if (con.send_gadgetron_configuration_file(std::string("isalive.xml")) != 0) {
		ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send XML configuration to the Gadgetron host")));
		return -1;
	}


	GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
			new GadgetContainerMessage<GadgetMessageIdentifier>();

	m1->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

	if (con.putq(m1) == -1) {
		ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to put CLOSE package on queue")));
		return -1;
	}

	con.wait();

	return 0;
}
