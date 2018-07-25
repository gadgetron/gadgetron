#pragma once

#include <mutex>
#include "gadgetron_home.h"
#include "gadgetron_matlab_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "engine.h"     // Matlab Engine header

//#include "ace/Synch.h"  // For the MatlabCommandServer
#include "ace/SOCK_Connector.h"
#include "ace/INET_Addr.h"


#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <boost/lexical_cast.hpp>
#include "mri_core_data.h"

// TODO:
//Create a debug option to use evalstring and get back the matlab output on every function call.
//Finish the image stuff
//Test on windows

extern std::mutex mutex_;

namespace Gadgetron{

class MatlabBufferGadget:
    public Gadget1<IsmrmrdReconData >
{
public:
    MatlabBufferGadget(): Gadget1<IsmrmrdReconData >()
    {
    }

    ~MatlabBufferGadget()
    {
    std::lock_guard<std::mutex> lock(mutex_);   
       // Close the Matlab engine
        GDEBUG("Closing down Matlab\n");
        engClose(engine_);
    }

    virtual int process(GadgetContainerMessage<IsmrmrdReconData> *);

protected:
    GADGET_PROPERTY(debug_mode, bool, "Debug mode", false);
    GADGET_PROPERTY(matlab_path, std::string, "Path to Matlab code", "");
    GADGET_PROPERTY(matlab_classname, std::string, "Name of Matlab gadget class", "");
    GADGET_PROPERTY(matlab_startcmd, std::string, "Matlab engine startup command", "matlab -nosplash");

    int process_config(ACE_Message_Block* mb)
    {
    	std::lock_guard<std::mutex> lock(mutex_);   
        std::string cmd;

        debug_mode_  = debug_mode.value();
        path_        = matlab_path.value();
        classname_   = matlab_classname.value();
	startcmd_    = matlab_startcmd.value();

        if (classname_.empty()) {
            GERROR("Missing Matlab Gadget classname in config!");
            return GADGET_FAIL;
        }

        GDEBUG("MATLAB Class Name : %s\n", classname_.c_str());


		// Open the Matlab Engine on the current host
		GDEBUG("Starting MATLAB engine with command: %s\n", startcmd_.c_str());
		if (!(engine_ = engOpen(startcmd_.c_str()))) {
			// TODO: error checking!
			GDEBUG("Can't start MATLAB engine\n");
		} else {
            // Prepare a buffer for collecting Matlab's output
            char matlab_buffer_[2049] = "\0";
            engOutputBuffer(engine_, matlab_buffer_, 2048);

	    // Add the necessary paths to the matlab environment
	    // Java matlab command server
	    std::string gadgetron_matlab_path = get_gadgetron_home().string() + "/share/gadgetron/matlab";
	    std::string add_path_cmd = std::string("addpath('") + gadgetron_matlab_path + std::string("');");
            // Gadgetron matlab scripts
	    engEvalString(engine_, add_path_cmd.c_str());
            // ISMRMRD matlab library
            engEvalString(engine_, "addpath(fullfile(getenv('ISMRMRD_HOME'), '/share/ismrmrd/matlab'));");
  
	    GDEBUG("%s", matlab_buffer_);
        }





        //char matlab_buffer_[2049] = "\0";
        char matlab_buffer_[20481] = "\0";
        engOutputBuffer(engine_, matlab_buffer_, 20480);

        // add user specified path for this gadget
        if (!path_.empty()) {
            cmd = "addpath('" + path_ + "');";
            send_matlab_command(cmd);
        }

        // Put the XML Header into the matlab workspace
        std::string xmlConfig = std::string(mb->rd_ptr());
        mxArray *xmlstring = mxCreateString(xmlConfig.c_str());
        engPutVariable(engine_, "xmlstring", xmlstring);

        // Instantiate the Matlab gadget object from the user specified class
        // Call matlab gadget's init method with the XML Header
        // and the user defined config method
        cmd = "matgadget = " + classname_ + "();";
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

            char matlab_buffer_[8193] = "\0";
            engOutputBuffer(engine_, matlab_buffer_, 8192);
            engEvalString(engine_, command.c_str());
	    if (debug_mode_) {
            GDEBUG("%s\n", matlab_buffer_);
	    }
            return GADGET_OK;

    }


    std::string path_;
    std::string classname_;
    std::string startcmd_;
    bool debug_mode_;

    Engine *engine_;
};

GADGET_FACTORY_DECLARE(MatlabBufferGadget);
}
