#include "ace/Log_Msg.h"
#include "ace/Get_Opt.h"
#include "ace/OS_NS_string.h"

#include "GadgetronConnector.h"
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "ImageWriter.h"
#include "HDF5ImageWriter.h"
#include "ImageAttribWriter.h"
#include "HDF5ImageAttribWriter.h"
#include "FileInfo.h"
#include "ismrmrd_dataset.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "BlobFileWriter.h"
#include "BlobFileWithAttribWriter.h"

#include <fstream>
#include <time.h>
#include <iomanip>

using namespace Gadgetron;
void print_usage()
{
    ACE_DEBUG((LM_INFO, ACE_TEXT("Usage: \n") ));
    ACE_DEBUG((LM_INFO, ACE_TEXT("mriclient -p <PORT>                                                                    (default 9002)\n") ));
    ACE_DEBUG((LM_INFO, ACE_TEXT("          -h <HOST>                                                                    (default localhost)\n") ));
    ACE_DEBUG((LM_INFO, ACE_TEXT("          -d <HDF5 DATA FILE>                                                          (default ./data.h5)\n") ));
    ACE_DEBUG((LM_INFO, ACE_TEXT("          -g <HDF5 DATA GROUP>                                                         (default /dataset)\n") ));
    ACE_DEBUG((LM_INFO, ACE_TEXT("          -c <GADGETRON CONFIG>                                                        (default default.xml)\n") ));
    ACE_DEBUG((LM_INFO, ACE_TEXT("          -l <LOOPS>                                                                   (default 1)\n") ));
    ACE_DEBUG((LM_INFO, ACE_TEXT("          -o <HDF5 OUT FILE>                                                           (out.h5)\n") ));
    ACE_DEBUG((LM_INFO, ACE_TEXT("          -G <HDF5 OUT GROUP>                                                          (default date and time)\n") ));
    ACE_DEBUG((LM_INFO, ACE_TEXT("          -F <OUT FILE FORMAT, 'h5 or hdf5' or 'hdf or analyze' >                      (default 'h5' format)\n") ));
}


std::string get_date_time_string()
{
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    std::stringstream str;
    str << timeinfo->tm_year+1900 << "-"
            << std::setw(2) << std::setfill('0') << timeinfo->tm_mon+1
            << "-"
            << std::setw(2) << std::setfill('0') << timeinfo->tm_mday
            << " "
            << std::setw(2) << std::setfill('0') << timeinfo->tm_hour
            << ":"
            << std::setw(2) << std::setfill('0') << timeinfo->tm_min
            << ":"
            << std::setw(2) << std::setfill('0') << timeinfo->tm_sec;

    std::string ret = str.str();

    return ret;
}

int ACE_TMAIN(int argc, ACE_TCHAR *argv[] )
{
    static const ACE_TCHAR options[] = ACE_TEXT(":p:h:d:x:c:l:o:g:G:F:");

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
    ACE_OS_String::strncpy(hdf5_out_file, "./out.h5", 1024);

    ACE_TCHAR hdf5_out_group[1024];

    std::string date_time = get_date_time_string();

    ACE_OS_String::strncpy(hdf5_out_group, date_time.c_str(), 1024);

    ACE_TCHAR out_format[128];
    ACE_OS_String::strncpy(out_format, "h5", 128);

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
        case 'o':
            ACE_OS_String::strncpy(hdf5_out_file, cmd_opts.opt_arg(), 1024);
            break;
        case 'G':
            ACE_OS_String::strncpy(hdf5_out_group, cmd_opts.opt_arg(), 1024);
            break;
        case 'F':
            ACE_OS_String::strncpy(out_format, cmd_opts.opt_arg(), 128);
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
    }

    boost::shared_ptr<ISMRMRD::Dataset> ismrmrd_dataset(new ISMRMRD::Dataset(hdf5_in_data_file,hdf5_in_group));
    std::string xml_config;
    ismrmrd_dataset->readHeader(xml_config);

    if (repetition_loops < 1) {
        ACE_DEBUG((LM_INFO, ACE_TEXT("Invalid number of repetition loops (%d).\n"), repetition_loops));
        print_usage();
        return -1;
    }

    ACE_DEBUG((LM_INFO, ACE_TEXT("  -- host            :      %s\n"), hostname));
    ACE_DEBUG((LM_INFO, ACE_TEXT("  -- port            :      %s\n"), port_no));
    ACE_DEBUG((LM_INFO, ACE_TEXT("  -- hdf5 file  in   :      %s\n"), hdf5_in_data_file));
    ACE_DEBUG((LM_INFO, ACE_TEXT("  -- hdf5 group in   :      %s\n"), hdf5_in_group));
    ACE_DEBUG((LM_INFO, ACE_TEXT("  -- conf            :      %s\n"), config_file));
    ACE_DEBUG((LM_INFO, ACE_TEXT("  -- loop            :      %d\n"), repetition_loops));
    ACE_DEBUG((LM_INFO, ACE_TEXT("  -- hdf5 file out   :      %s\n"), hdf5_out_file));
    ACE_DEBUG((LM_INFO, ACE_TEXT("  -- hdf5 group out  :      %s\n"), hdf5_out_group));
    ACE_DEBUG((LM_INFO, ACE_TEXT("  -- out format      :      %s\n"), out_format));

    std::string prefix;
    std::string out_format_str(out_format);

    for (int i = 0; i < repetition_loops; i++)
    {
        if ( repetition_loops > 1 )
        {
            std::ostringstream ostr;
            ostr << "MriClient_Run" << i;
            prefix = ostr.str();
        }

        GadgetronConnector con;

        //con.register_writer(GADGET_MESSAGE_ACQUISITION, new GadgetAcquisitionMessageWriter());
        con.register_writer(GADGET_MESSAGE_ISMRMRD_ACQUISITION, new GadgetIsmrmrdAcquisitionMessageWriter());
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT, new HDF5ImageWriter<ACE_UINT16>(std::string(hdf5_out_file), std::string(hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT, new HDF5ImageWriter<float>(std::string(hdf5_out_file), std::string(hdf5_out_group)));
        con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT, new HDF5ImageWriter< std::complex<float> >(std::string(hdf5_out_file), std::string(hdf5_out_group)));

        if ( (out_format_str == "analyze") || (out_format_str == "hdr") )
        {
            con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT, new AnalyzeImageAttribWriter<ACE_UINT16>(prefix));
            con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT, new AnalyzeImageAttribWriter<float>(prefix));
            con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT, new AnalyzeComplexImageAttribWriter< std::complex<float> >(prefix));
        }
        else
        {
            con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT, new HDF5ImageAttribWriter<ACE_UINT16>(std::string(hdf5_out_file), std::string(hdf5_out_group), prefix));
            con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT, new HDF5ImageAttribWriter<float>(std::string(hdf5_out_file), std::string(hdf5_out_group), prefix));
            con.register_reader(GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT, new HDF5ImageAttribWriter< std::complex<float> >(std::string(hdf5_out_file), std::string(hdf5_out_group), prefix));
        }

        con.register_reader(GADGET_MESSAGE_DICOM, new BlobFileWriter(std::string(hdf5_out_group), std::string("dcm")));
        con.register_reader(GADGET_MESSAGE_DICOM_WITHNAME, new BlobFileWithAttribWriter(std::string(), std::string("dcm")));

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

        if (con.send_gadgetron_parameters(xml_config) != 0) {
            ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send XML parameters to the Gadgetron host")));
            return -1;
        }

        unsigned long acquisitions = ismrmrd_dataset->getNumberOfAcquisitions();//HDF5GetLengthOfFirstDimension(hdf5_in_data_file, hdf5_data_varname.c_str());

        for (unsigned long int i = 0; i < acquisitions; i++) {
            GadgetContainerMessage<ISMRMRD::Acquisition>* acq = new GadgetContainerMessage<ISMRMRD::Acquisition>();
	    ismrmrd_dataset->readAcquisition(i, *acq->getObjectPtr());
            
            GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
                    new GadgetContainerMessage<GadgetMessageIdentifier>();

            m1->getObjectPtr()->id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;

            m1->cont(acq);

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
