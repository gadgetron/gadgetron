#include "GadgetronConnector.h"
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "DependencyQueryReader.h"

#include <ace/Log_Msg.h>
#include <ace/Get_Opt.h>
#include <ace/OS_NS_string.h>

#include <fstream>
#include <time.h>
#include <iomanip>
#include <sstream>

using namespace Gadgetron;

static void usage()
{
    using namespace std;
    std::ostringstream outs;

    outs << "Query the gadgetron server for the stored dependency measurements" << endl;
    outs << "gtdependencyquery   -p <PORT>                      (default 9002)" << endl;
    outs << "                    -h <HOST>                      (default localhost)" << endl;
    outs << "                    -o <Query out file>            (default dependency.xml)" << endl;
    outs << std::ends; 

    GDEBUG_STREAM(outs.str());
}

int ACE_TMAIN(int argc, ACE_TCHAR *argv[] )
{
    GadgetronConnector con;

    std::string host("localhost");
    std::string port("9002");
    std::string out("dependency.xml");

    static const ACE_TCHAR options[] = ACE_TEXT(":p:h:o:");

    ACE_Get_Opt cmd_opts(argc, argv, options);

    int option;
    while ((option = cmd_opts()) != EOF)
    {
        switch (option) {
        case 'p':
            port = std::string(cmd_opts.opt_arg());
            break;
        case 'h':
            host = std::string(cmd_opts.opt_arg());
            break;
        case 'o':
            out = std::string(cmd_opts.opt_arg());
            break;
        case ':':
            usage();
            GERROR("-%c requires an argument.\n", cmd_opts.opt_opt());
	    return -1;
            break;
        default:
            usage();
            GERROR("Command line parse error\n");
	    return -1;
            break;
        }
    }

    if (con.open(host,port) != 0)
    {
      GERROR("Unable to connect to the Gadgetron host\n");
      return -1;
    }

    // need to register a reader
    con.register_reader(GADGET_MESSAGE_DEPENDENCY_QUERY, new DependencyQueryReader(out));

    //Tell Gadgetron which XML configuration to run.
    if (con.send_gadgetron_configuration_file(std::string("gtquery.xml")) != 0)
    {
      GERROR("Unable to send XML configuration to the Gadgetron host\n");
      return -1;
    }

    GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
            new GadgetContainerMessage<GadgetMessageIdentifier>();

    m1->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

    if (con.putq(m1) == -1)
    {
      GERROR("Unable to put CLOSE package on queue\n");
      return -1;
    }

    con.wait();

    return 0;
}
