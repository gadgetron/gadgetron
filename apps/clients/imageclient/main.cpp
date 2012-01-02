#include <QApplication>
#include <fstream>

#include "GadgetSocketSender.h" //We need this for now for the GadgetAcquisitionWriter
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "../mriclient/ImageWriter.h"
#include "FileInfo.h"

#include "mainWindow.h"

void print_usage()
{
	ACE_DEBUG((LM_INFO, ACE_TEXT("Usage: \n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("imageclient -p <PORT>                      (default 9002)       \n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("            -h <HOST>                      (default localhost)  \n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("            -c <GADGETRON CONFIG>                               \n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("            -l <LOOPS>                     (default 1)          \n") ));
}

int main(int argc, char *argv[])
{
	// Parse command line parameters
	//
	//

	static const ACE_TCHAR options[] = ACE_TEXT(":p:h:c:l:");

	ACE_Get_Opt cmd_opts(argc, argv, options);

	ACE_TCHAR port_no[1024];
	ACE_OS_String::strncpy(port_no, "9002", 1024);

	ACE_TCHAR hostname[1024];
	ACE_OS_String::strncpy(hostname, "localhost", 1024);

	ACE_TCHAR config_file[1024];
	ACE_OS_String::strncpy(config_file, "", 1024);

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

	if (repetition_loops < 1) {
		ACE_DEBUG((LM_INFO, ACE_TEXT("Invalid number of repetition loops (%d).\n"), repetition_loops));
		print_usage();
		return -1;
	}

	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- host:            %s\n"), hostname));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- port:            %s\n"), port_no));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- conf:            %s\n"), config_file));
	ACE_DEBUG((LM_INFO, ACE_TEXT("  -- loop:            %d\n"), repetition_loops));

	ACE_DEBUG(( LM_INFO, ACE_TEXT("Gadgetron image sender started\n") ));

	// Open QT interface
	//
	QApplication app(argc, argv);
	MainWindow mw;
	mw.show();

	// Pass command line settings
	//
	mw.set_hostname(hostname);
	mw.set_port_no(port_no);
	mw.set_config(config_file);
	mw.set_num_repetitions(repetition_loops);

	// Establish connection to the Gadgetron
	if( mw.setup_gadgetron_connection() < 0) return 1;

	// Everything else is triggered from the GUI
	//
	return app.exec();
}
