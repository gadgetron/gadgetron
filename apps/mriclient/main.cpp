#include "ace/Log_Msg.h"
#include "ace/Get_Opt.h"
#include "ace/OS_NS_string.h"

#include "GadgetronConnector.h"
#include "GadgetSocketSender.h" //We need this for now for the GadgetAcquisitionWriter
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "ImageWriter.h"
#include "FileInfo.h"

#include <fstream>

void print_usage()
{
	ACE_DEBUG((LM_INFO, ACE_TEXT("Usage: \n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("mriclient -p <PORT>                      (default 9002)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -h <HOST>                      (default localhost)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -d <DATA FILE>                 (default ./data.dat)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -x <PARAMETER FILE (XML)>      (default ./parameters.xml)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -c <GADGETRON CONFIG>          (default default.xml)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -l <LOOPS>                     (default 1)\n") ));
}

int ACE_TMAIN(int argc, ACE_TCHAR *argv[] )
{
	static const ACE_TCHAR options[] = ACE_TEXT(":p:h:d:x:c:l:");

	ACE_Get_Opt cmd_opts(argc, argv, options);

	ACE_TCHAR port_no[1024];
	ACE_OS_String::strncpy(port_no, "9002", 1024);

	ACE_TCHAR hostname[1024];
	ACE_OS_String::strncpy(hostname, "localhost", 1024);

	ACE_TCHAR data_file[4096];
	ACE_OS_String::strncpy(data_file, "./data.dat", 4096);

	ACE_TCHAR parameter_file[4096];
	ACE_OS_String::strncpy(parameter_file, "./parameters.xml", 4096);

	ACE_TCHAR config_file[1024];
	ACE_OS_String::strncpy(config_file, "default.xml", 1024);

	int repetition_loops = 1;

	int option;
	while ((option = cmd_opts()) != EOF) {
		switch (option) {
		case 'p':
			ACE_OS_String::strncpy(port_no, cmd_opts.opt_arg(), 1024);
			break;
		case 'h':
			ACE_OS_String::strncpy(hostname, cmd_opts.opt_arg(), 1024);
			break;
		case 'd':
			ACE_OS_String::strncpy(data_file, cmd_opts.opt_arg(), 4096);
			break;
		case 'x':
			ACE_OS_String::strncpy(parameter_file, cmd_opts.opt_arg(), 4096);
			break;
		case 'c':
			ACE_OS_String::strncpy(config_file, cmd_opts.opt_arg(), 1024);
			break;
		case 'l':
			repetition_loops = ACE_OS::atoi(cmd_opts.opt_arg());
			break;
		case ':':
			print_usage();
			ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("-%c requires an argument.\n"), cmd_opts.opt_opt()),-1);
			break;
		default:
			print_usage();
			ACE_ERROR_RETURN( (LM_ERROR, ACE_TEXT("Command line parse error\n")), -1);
			break;
		}
	}

	ACE_DEBUG(( LM_INFO, ACE_TEXT("Gadgetron MRI Data Sender\n") ));

	//Let's check if the files exist:
	if (!FileInfo(std::string(data_file)).exists()) {
		ACE_DEBUG((LM_INFO, ACE_TEXT("Data file %s does not exist.\n"), data_file));
		print_usage();
		return -1;
	}

	if (!FileInfo(std::string(parameter_file)).exists()) {
		ACE_DEBUG((LM_INFO, ACE_TEXT("Parameter file %s does not exist.\n"), parameter_file));
		print_usage();
		return -1;
	}

	if (repetition_loops < 1) {
		ACE_DEBUG((LM_INFO, ACE_TEXT("Invalid number of repetition loops (%d).\n"), repetition_loops));
		print_usage();
		return -1;
	}

	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- host:            %s\n"), hostname));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- port:            %s\n"), port_no));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- data:            %s\n"), data_file));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- parm:            %s\n"), parameter_file));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- conf:            %s\n"), config_file));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- loop:            %d\n"), repetition_loops));

	for (int i = 0; i < repetition_loops; i++) {

		GadgetronConnector con;

		con.register_writer(GADGET_MESSAGE_ACQUISITION, new GadgetAcquisitionMessageWriter());
		con.register_reader(GADGET_MESSAGE_IMAGE_REAL_USHORT, new ImageWriter<ACE_UINT16>());
		con.register_reader(GADGET_MESSAGE_IMAGE_REAL_FLOAT, new ImageWriter<float>());
		con.register_reader(GADGET_MESSAGE_IMAGE_CPLX_FLOAT, new ImageWriter< std::complex<float> >());

		//Open a connection with the gadgetron
		if (con.open(std::string(hostname),std::string(port_no)) != 0) {
			ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to connect to the Gadgetron host")));
			return -1;
		}

		//Tell Gadgetron which XML configuration to run.
		if (con.send_gadgetron_configuration(std::string(config_file)) != 0) {
			ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send XML configuration to the Gadgetron host")));
			return -1;
		}

		//Read XML parameter file from disk
		std::ifstream ifs(parameter_file);
		std::string xmlstring = std::string(std::istreambuf_iterator<char>(ifs),std::istreambuf_iterator<char>());

		if (con.send_gadgetron_parameters(xmlstring) != 0) {
			ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send XML parameters to the Gadgetron host")));
			return -1;
		}

		//We are now ready to send data...
		std::ifstream is;
		is.open (data_file, ios::binary );
		is.seekg (0, ios::end);
		size_t length = is.tellg();
		is.seekg (0, ios::beg);

		while ((length-is.tellg()) > sizeof(GadgetMessageAcquisition)) {
			GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
					new GadgetContainerMessage<GadgetMessageIdentifier>();

			m1->getObjectPtr()->id = GADGET_MESSAGE_ACQUISITION;

			GadgetContainerMessage<GadgetMessageAcquisition>* m2 =
					new GadgetContainerMessage<GadgetMessageAcquisition>();

			//Read acquisition header from disk
			is.read(reinterpret_cast<char*>(m2->getObjectPtr()), sizeof(GadgetMessageAcquisition));

			std::vector<unsigned int> dimensions(2);
			dimensions[0] = m2->getObjectPtr()->samples;
			dimensions[1] = m2->getObjectPtr()->channels;

			GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 =
					new GadgetContainerMessage< hoNDArray< std::complex< float> > >();

			if (!m3->getObjectPtr()->create(&dimensions)) {
				ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to create storage for NDArray.\n")));
				ACE_DEBUG((LM_ERROR, ACE_TEXT("Requested dimensions were (%d,%d)"), dimensions[0], dimensions[1]));
				return -1;
			}

			//Read data from disk
			is.read(reinterpret_cast<char*>(m3->getObjectPtr()->get_data_ptr()),sizeof(float)*2*m3->getObjectPtr()->get_number_of_elements());

			//Chain the message block together.
			m1->cont(m2);
			m2->cont(m3);

			if (con.putq(m1) == -1) {
				ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to put data package on queue")));
				return -1;
			}
		}

		is.close();

		GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
				new GadgetContainerMessage<GadgetMessageIdentifier>();

		m1->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

		if (con.putq(m1) == -1) {
			ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to put CLOSE package on queue")));
			return -1;
		}

		con.wait();
	}

	return 0;
}
