#pragma once

#include "gadgetron_matlab_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "engine.h"     // Matlab Engine header

#include "ace/Synch.h"  // For the MatlabCommandServer
#include "ace/SOCK_Connector.h"
#include "ace/INET_Addr.h"


#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <boost/lexical_cast.hpp>
#include "mri_core_data.h"
// TODO:
//Make the port option work so that we can have multiple matlabs running, each with its own command server.
//Create a debug option to use evalstring and get back the matlab output on every function call.
//Finish the image stuff
//Is there a better way to kill the command server?
//Test on windows


namespace Gadgetron{

class MatlabBufferGadget:
    public Gadget1<IsmrmrdReconData >
{
public:
    MatlabBufferGadget(): Gadget1<IsmrmrdReconData >()
    {
        // Open the Matlab Engine on the current host
        GDEBUG("Starting MATLAB engine\n");
        if (!(engine_ = engOpen("matlab -nosplash -nodesktop"))) {
            // TODO: error checking!
            GDEBUG("Can't start MATLAB engine\n");
        } else {
            // Add ISMRMRD Java bindings jar to Matlab's path
            // TODO: this should be in user's Matlab path NOT HERE

            // Prepare a buffer for collecting Matlab's output
            char matlab_buffer_[2049] = "\0";
            engOutputBuffer(engine_, matlab_buffer_, 2048);

	    // Add the necessary paths to the matlab environment
	    // Java matlab command server
            engEvalString(engine_, "javaaddpath(fullfile(getenv('GADGETRON_HOME'), 'matlab'));");
            // Gadgetron matlab scripts
            engEvalString(engine_, "addpath(fullfile(getenv('GADGETRON_HOME'), 'matlab'));");
            // ISMRMRD matlab library
            engEvalString(engine_, "addpath(fullfile(getenv('ISMRMRD_HOME'), 'matlab'));");


	    GDEBUG("%s", matlab_buffer_);
        }
    }

    ~MatlabBufferGadget()
    {
        char matlab_buffer_[2049] = "\0";
        engOutputBuffer(engine_, matlab_buffer_, 2048);
	// Stop the Java Command server
        // send the stop signal to the command server and
        //  wait a bit for it to shut down cleanly.
        GDEBUG("Closing down the Matlab Command Server\n");
	engEvalString(engine_, "M.notifyEnd(); pause(1);");
        engEvalString(engine_, "clear java;");
        GDEBUG("%s", matlab_buffer_);
        // Close the Matlab engine
        GDEBUG("Closing down Matlab\n");
        engClose(engine_);
    }

    virtual int process(GadgetContainerMessage<IsmrmrdReconData> *);

protected:

    int process_config(ACE_Message_Block* mb)
    {
        std::string cmd;

        debug_mode_  = this->get_int_value("debug_mode");
        path_        = this->get_string_value("matlab_path");
        classname_   = this->get_string_value("matlab_classname");
        command_server_port_ = this->get_int_value("matlab_port");

        GDEBUG("MATLAB Class Name : %s\n", classname_.get()->c_str());

        //char matlab_buffer_[2049] = "\0";
        char matlab_buffer_[20481] = "\0";
        engOutputBuffer(engine_, matlab_buffer_, 20480);

   	// Instantiate the Java Command server
        // TODO: we HAVE to pause in Matlab to allow the java command server thread to start
        cmd = "M = MatlabCommandServer(" + boost::lexical_cast<std::string>(command_server_port_) +
                "); M.start(); pause(1);";
	engEvalString(engine_, cmd.c_str());
        GDEBUG("%s", matlab_buffer_);

        // add user specified path for this gadget
        if (!path_->empty()) {
            cmd = "addpath('" + *path_ + "');";
            send_matlab_command(cmd);
        }

        // Put the XML Header into the matlab workspace
        std::string xmlConfig = std::string(mb->rd_ptr());
        mxArray *xmlstring = mxCreateString(xmlConfig.c_str());
        engPutVariable(engine_, "xmlstring", xmlstring);

        // Instantiate the Matlab gadget object from the user specified class
        // Call matlab gadget's init method with the XML Header
        // and the user defined config method
        cmd = "matgadget = " + *classname_ + "();";
        cmd += "matgadget.init(xmlstring); matgadget.config();";
        if (send_matlab_command(cmd) != GADGET_OK) {
            GDEBUG("Failed to send matlab command.\n");
            return GADGET_FAIL;
        }

	mxDestroyArray(xmlstring);

        return GADGET_OK;
    }

    int send_matlab_command(std::string& command)
    {

        if (debug_mode_) {
            char matlab_buffer_[2049] = "\0";
            engOutputBuffer(engine_, matlab_buffer_, 2048);
            engEvalString(engine_, command.c_str());
            GDEBUG("%s\n", matlab_buffer_);
            return GADGET_OK;
        }
        else {
            ACE_SOCK_Stream client_stream;
            ACE_INET_Addr remote_addr(command_server_port_, "localhost");
            ACE_SOCK_Connector connector;

            if (connector.connect(client_stream, remote_addr) == -1) {
                GDEBUG("Connection failed\n");
                return GADGET_FAIL;
            }

            ACE_Time_Value timeout(10);
            if (client_stream.send_n(command.c_str(), command.size(), &timeout) == -1) {
                GDEBUG("Error in send_n\n");
                client_stream.close();
                return GADGET_FAIL;
            }

            if (client_stream.close () == -1){
                GDEBUG("Error in close\n");
                return GADGET_FAIL;
            }
            return GADGET_OK;
        }
    }


    boost::shared_ptr<std::string> path_;
    boost::shared_ptr<std::string> classname_;
    int command_server_port_;
    int debug_mode_;

    Engine *engine_;
};

GADGET_FACTORY_DECLARE(MatlabBufferGadget);
}
