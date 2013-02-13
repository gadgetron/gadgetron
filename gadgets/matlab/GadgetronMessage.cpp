//#include "ace/FIFO_Send_Msg.h"
#include "ace/OS_NS_stropts.h"
#include "ace/Log_Msg.h"
#include "ace/OS_NS_stdio.h"
#include "ace/OS_NS_unistd.h"
#include "ace/SOCK_Connector.h"

#include "MatlabCommunicator.h"

#include <stdio.h>
#include "mex.h"


void  mexFunction(int nlhs,mxArray *[],int nrhs,const mxArray *prhs[])
{
	char *str;

	/* Check for proper number of input and output arguments */
	if (nrhs != 1) {
		mexErrMsgIdAndTxt( "MATLAB:mexatexit:invalidNumInputs",
				"One input argument required.");
	}
	if (nlhs > 1){
		mexErrMsgIdAndTxt( "MATLAB:mexatexit:maxlhs",
				"Too many output arguments.");
	}

	/* Check to be sure input is of type char */
	if (!(mxIsChar(prhs[0]))){
		mexErrMsgIdAndTxt( "MATLAB:mexatexit:inputNotString",
				"Input must be of type string.\n.");
	}

	str = mxArrayToString(prhs[0]);
	size_t len = mxGetNumberOfElements(prhs[0]);

	mexPrintf("Sending string: %s [%d].\n", str, len);


	mexPrintf("Allocating Buffer\n");
	ACE_Str_Buf msg (str, len);

	mexPrintf("Opening FIFO\n");

	ACE_SOCK_Stream cli_stream;
	ACE_SOCK_Connector con;
	int i;
	const ACE_INET_Addr addr("9003");
	if (con.connect (cli_stream, addr) == -1) {
		mexErrMsgTxt("Error Opening PIPE to Gadgetron");
	}


	if (cli_stream.send (msg.buf, msg.len) != msg.len) {
		mexErrMsgTxt("Error Sending to Gadgetron");
	}

	if (cli_stream.close () == -1){
		mexErrMsgTxt("Error shutting down stream");
	}


	mexPrintf("Done sending data\n");

	mxFree(str);
	return;
}
