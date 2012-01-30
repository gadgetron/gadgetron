#include "ace/Log_Msg.h"
#include "ace/Get_Opt.h"
#include "ace/OS_NS_string.h"

#include "GadgetronConnector.h"
#include "GadgetSocketSender.h" //We need this for now for the GadgetAcquisitionWriter
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "ImageWriter.h"
#include "HDF5ImageWriter.h"
#include "FileInfo.h"
#include "mri_hdf5_io.h"

#include <fstream>
#include <time.h>

void print_usage()
{
	ACE_DEBUG((LM_INFO, ACE_TEXT("Usage: \n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("mriclient -p <PORT>                      (default 9002)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -h <HOST>                      (default localhost)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -d <HDF5 DATA FILE>            (default ./data.h5)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -g <HDF5 Group>                (default /dataset)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -c <GADGETRON CONFIG>          (default default.xml)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -l <LOOPS>                     (default 1)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -5 <HDF5 out file name>        (no default)\n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("          -G <HDF5 out group name>       (default date and time)\n") ));
}

int ACE_TMAIN(int argc, ACE_TCHAR *argv[] )
{
	static const ACE_TCHAR options[] = ACE_TEXT(":p:h:d:x:c:l:5:g:G:");

	ACE_Get_Opt cmd_opts(argc, argv, options);

	ACE_TCHAR port_no[1024];
	ACE_OS_String::strncpy(port_no, "9002", 1024);

	ACE_TCHAR hostname[1024];
	ACE_OS_String::strncpy(hostname, "localhost", 1024);

	ACE_TCHAR hdf5_in_data_file[4096];
	ACE_OS_String::strncpy(hdf5_in_data_file, "./data.h5", 4096);

	ACE_TCHAR hdf5_in_group[4096];
	ACE_OS_String::strncpy(hdf5_in_group, "/dataset", 4096);

	ACE_TCHAR config_file[1024];
	ACE_OS_String::strncpy(config_file, "default.xml", 1024);

	bool save_hdf5 = false;

	ACE_TCHAR hdf5_out_file[1024];
	ACE_TCHAR hdf5_out_group[1024];

	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	std::stringstream str;
	str << timeinfo->tm_year+1900 << "-" << timeinfo->tm_mon+1 << "-" << timeinfo->tm_mday
		<< " " << timeinfo->tm_hour << ":" << timeinfo->tm_min << ":" << timeinfo->tm_sec;

	ACE_OS_String::strncpy(hdf5_out_group, str.str().c_str(), 1024);

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
			ACE_OS_String::strncpy(hdf5_in_data_file, cmd_opts.opt_arg(), 4096);
			break;
		case 'g':
			ACE_OS_String::strncpy(hdf5_in_group, cmd_opts.opt_arg(), 4096);
			break;
		case 'c':
			ACE_OS_String::strncpy(config_file, cmd_opts.opt_arg(), 1024);
			break;
		case 'l':
			repetition_loops = ACE_OS::atoi(cmd_opts.opt_arg());
			break;
		case '5':
			ACE_OS_String::strncpy(hdf5_out_file, cmd_opts.opt_arg(), 1024);
			save_hdf5 = true;
			break;
		case 'G':
			ACE_OS_String::strncpy(hdf5_out_group, cmd_opts.opt_arg(), 1024);
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
	std::string hdf5_xml_varname = std::string(hdf5_in_group) + std::string("/xml");
	std::string hdf5_data_varname = std::string(hdf5_in_group) + std::string("/data");

	if (!FileInfo(std::string(hdf5_in_data_file)).exists()) {
		ACE_DEBUG((LM_INFO, ACE_TEXT("Data file %s does not exist.\n"), hdf5_in_data_file));
		print_usage();
		return -1;
	} else {
		boost::shared_ptr<H5File> f = OpenHF5File(hdf5_in_data_file);
		if (!HDF5LinkExists(f.get(), hdf5_xml_varname.c_str())) {
			  ACE_DEBUG((LM_INFO, ACE_TEXT("HDF5 variable \"%s\" does not exists.\n"), hdf5_xml_varname.c_str()));
			  print_usage();
			  return -1;
		}
		if (!HDF5LinkExists(f.get(), hdf5_data_varname.c_str())) {
			  ACE_DEBUG((LM_INFO, ACE_TEXT("HDF5 variable \"%s\" does not exists.\n"), hdf5_data_varname.c_str()));
			  print_usage();
			  return -1;
		}
	}


	if (repetition_loops < 1) {
		ACE_DEBUG((LM_INFO, ACE_TEXT("Invalid number of repetition loops (%d).\n"), repetition_loops));
		print_usage();
		return -1;
	}

	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- host:            %s\n"), hostname));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- port:            %s\n"), port_no));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- data (in):       %s\n"), hdf5_in_data_file));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- group (in):      %s\n"), hdf5_in_group));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- conf:            %s\n"), config_file));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- loop:            %d\n"), repetition_loops));

	if (save_hdf5) {
		ACE_DEBUG((LM_INFO, ACE_TEXT("  -- hdf5 (file)    :      %s\n"), hdf5_out_file));
		ACE_DEBUG((LM_INFO, ACE_TEXT("  -- hdf5 (group)   :      %s\n"), hdf5_out_group));
	}

	for (int i = 0; i < repetition_loops; i++) {

		GadgetronConnector con;

		con.register_writer(GADGET_MESSAGE_ACQUISITION, new GadgetAcquisitionMessageWriter());
		if (save_hdf5) {
			con.register_reader(GADGET_MESSAGE_IMAGE_REAL_USHORT, new HDF5ImageWriter<ACE_UINT16>(std::string(hdf5_out_file), std::string(hdf5_out_group)));
			con.register_reader(GADGET_MESSAGE_IMAGE_REAL_FLOAT, new HDF5ImageWriter<float>(std::string(hdf5_out_file), std::string(hdf5_out_group)));
			con.register_reader(GADGET_MESSAGE_IMAGE_CPLX_FLOAT, new HDF5ImageWriter< std::complex<float> >(std::string(hdf5_out_file), std::string(hdf5_out_group)));

		} else {
			con.register_reader(GADGET_MESSAGE_IMAGE_REAL_USHORT, new ImageWriter<ACE_UINT16>());
			con.register_reader(GADGET_MESSAGE_IMAGE_REAL_FLOAT, new ImageWriter<float>());
			con.register_reader(GADGET_MESSAGE_IMAGE_CPLX_FLOAT, new ImageWriter< std::complex<float> >());
		}

		//Open a connection with the gadgetron
		if (con.open(std::string(hostname),std::string(port_no)) != 0) {
			ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to connect to the Gadgetron host")));
			return -1;
		}

		//Tell Gadgetron which XML configuration to run.
		if (con.send_gadgetron_configuration_file(std::string(config_file)) != 0) {
			ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send XML configuration to the Gadgetron host")));
			return -1;
		}

		//Read XML parameter file from disk
		boost::shared_ptr< hoNDArray<char> > xml_array = hdf5_read_array_slice<char>(hdf5_in_data_file, hdf5_xml_varname.c_str());

		std::string xmlstring(xml_array->get_data_ptr());

		if (con.send_gadgetron_parameters(xmlstring) != 0) {
			ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send XML parameters to the Gadgetron host")));
			return -1;
		}

		unsigned long acquisitions = HDF5GetLengthOfFirstDimension(hdf5_in_data_file, hdf5_data_varname.c_str());

		for (unsigned long int i = 0; i < acquisitions; i++) {

			header_data_struct<GadgetMessageAcquisition, std::complex<float> > acquisition;

			{

				HDF5Exclusive lock; //This will ensure threadsafe access to HDF5
				acquisition = hdf5_read_struct_with_data< GadgetMessageAcquisition, std::complex<float> >(hdf5_in_data_file, hdf5_data_varname.c_str(), i);
			}


			GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
					new GadgetContainerMessage<GadgetMessageIdentifier>();

			m1->getObjectPtr()->id = GADGET_MESSAGE_ACQUISITION;

			GadgetContainerMessage<GadgetMessageAcquisition>* m2 =
					new GadgetContainerMessage<GadgetMessageAcquisition>();

			//We have to copy this data, because the header_data_struct contains boost::shared_ptr and those objects will get destroyed.
			//There is no clean way to take the ownership from a boost::shared_ptr and hand it to an ACE_Message_Block
			memcpy(m2->getObjectPtr(), acquisition.h.get(), sizeof(GadgetMessageAcquisition));

			std::vector<unsigned int> dimensions(2);
			dimensions[0] = m2->getObjectPtr()->samples;
			dimensions[1] = m2->getObjectPtr()->channels;

			GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 =
					new GadgetContainerMessage< hoNDArray< std::complex< float> > >();

			//We have to copy, see above.
			*m3->getObjectPtr() = *acquisition.d.get();

			//Chain the message block together.
			m1->cont(m2);
			m2->cont(m3);

			if (con.putq(m1) == -1) {
				ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to put data package on queue")));
				return -1;
			}
		}

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
