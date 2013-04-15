#pragma once

#include "gadgetronmatlab_export.h"
#include "Gadget.h"
#include "Gadgetron.h"
#include "hoNDArray.h"
#include "ismrmrd.h"

#include "engine.h"     // Matlab Engine header

#include "ace/Synch.h"  // For the MatlabCommandServer
#include "ace/SOCK_Connector.h"
#include "ace/INET_Addr.h"


#include <stdio.h>
#include <stdlib.h>

#include <complex>
#include <boost/lexical_cast.hpp>

// TODO:
//Make the port option work so that we can have multiple matlabs running, each with its own command server.
//Create a debug option to use evalstring and get back the matlab output on every function call.
//Finish the image stuff
//Is there a better way to kill the command server?
//Test on windows


template <class T> class MatlabGadget :
    public Gadget2<T, hoNDArray< std::complex<float> > >
{
public:
    MatlabGadget(): Gadget2<T, hoNDArray< std::complex<float> > >()
    {
        // Open the Matlab Engine on the current host
        GADGET_DEBUG1("Starting MATLAB engine\n");
        if (!(engine_ = engOpen("matlab -nosplash -nodesktop"))) {
            // TODO: error checking!
            GADGET_DEBUG1("Can't start MATLAB engine\n");
        } else {
            // Add ISMRMRD Java bindings jar to Matlab's path
            // TODO: this should be in user's Matlab path NOT HERE
	  
            // Prepare a buffer for collecting Matlab's output
            char buffer[2049] = "\0";
            engOutputBuffer(engine_, buffer, 2048);
            engEvalString(engine_, "javaaddpath(fullfile(getenv('ISMRMRD_HOME'),'lib','ismrmrd.jar'));");
            engEvalString(engine_, "javaaddpath(fullfile(getenv('GADGETRON_HOME'), 'matlab'));");
            engEvalString(engine_, "addpath(fullfile(getenv('GADGETRON_HOME'), 'matlab'));");
            engEvalString(engine_, "addpath(fullfile(getenv('ISMRMRD_HOME'), 'matlab'));");

            // Load the Java JNI library
            engEvalString(engine_, "org.ismrm.ismrmrd.JNILibLoader.load();");

	    GADGET_DEBUG2("%s", buffer);
        }
    }

    ~MatlabGadget()
    {
        char buffer[2049] = "\0";
        engOutputBuffer(engine_, buffer, 2048);
	// Stop the Java Command server
        // send the stop signal to the command server and
        //  wait a bit for it to shut down cleanly.
        GADGET_DEBUG1("Closing down the Matlab Command Server\n");
	engEvalString(engine_, "M.notifyEnd(); pause(2);");
        engEvalString(engine_, "clear java;");
        GADGET_DEBUG2("%s", buffer);
        // Close the Matlab engine
        GADGET_DEBUG1("Closing down Matlab\n");
        engClose(engine_);
    }

protected:

    int process_config(ACE_Message_Block* mb)
    {

        std::string cmd = "";

        path_        = this->get_string_value("path");
        classname_   = this->get_string_value("matlab_classname");
        //mport        = this->get_string_value("matlab_port");

        GADGET_DEBUG2("MATLAB Path    : %s\n", path_.get()->c_str());
        GADGET_DEBUG2("MATLAB Class Name : %s\n", classname_.get()->c_str());

        // Set up buffer for catching Matlab output
        char buffer[2049] = "\0";
        engOutputBuffer(engine_, buffer, 2048);

   	// Instantiate the Java Command server
        // TODO: make the command server run on a user specified port
	engEvalString(engine_, "M = MatlabCommandServer(); M.start();");
        //GADGET_DEBUG2("%s", buffer);

        // add user specified path for this gadget
        if (path_->length() > 0) {
            cmd = "addpath(" + *path_.get() + ");";
            //engEvalString(engine_, cmd.c_str());
            //GADGET_DEBUG2("%s", buffer);
            send_matlab_command(cmd);
        }

        // Put the XML Header into the matlab workspace
        std::string xmlConfig = std::string(mb->rd_ptr());
        mxArray *xmlstring = mxCreateString(xmlConfig.c_str());
        engPutVariable(engine_, "xmlstring", xmlstring);

        // Instantiate the Matlab gadget object from the user specified class
        // Call matlab gadget's init method with the XML Header
        // and the user definined config method
        cmd = "matgadget = " + *classname_.get() + "();";
        cmd += "matgadget.init(xmlstring); matgadget.config();";
        //engEvalString(engine_, cmd.c_str());
        //GADGET_DEBUG2("%s", buffer);
        send_matlab_command(cmd);

	mxDestroyArray(xmlstring);

        return GADGET_OK;
    }

    int send_matlab_command(std::string& command)
    {
      ACE_SOCK_Stream client_stream_;
      ACE_INET_Addr remote_addr_("3000","localhost");
      ACE_SOCK_Connector connector_;

      if (connector_.connect (client_stream_, remote_addr_) == -1) {
	std::cout << "Connection failed" << endl;
        return -1;
      }
      else  {
        //std::cout<<"Connected to " << remote_addr_.get_host_name() << ":" << remote_addr_.get_port_number() << endl;
      }

      ACE_Time_Value timeout(10);
      if (client_stream_.send_n(command.c_str(), command.size(), &timeout) == -1) {
	std::cout << "Error in send_n" << std::endl;
      }
      else {
	//std::cout << "Sent " << command << std::endl;
      }

      if (client_stream_.close () == -1){
	std::cout << "Error in close" << std::endl;
	return -1;
      }
      return 0;
    }

    boost::shared_ptr<std::string> path_;
    boost::shared_ptr<std::string> classname_;

    Engine *engine_;

};



class EXPORTGADGETSMATLAB AcquisitionMatlabGadget :
    public MatlabGadget<ISMRMRD::AcquisitionHeader>
{
    public:
        GADGET_DECLARE(AcquisitionMatlabGadget);

        int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
                GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

};

class EXPORTGADGETSMATLAB ImageMatlabGadget :
    public MatlabGadget<ISMRMRD::ImageHeader>
{
    public:
        GADGET_DECLARE(ImageMatlabGadget);

        int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
                GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

};
