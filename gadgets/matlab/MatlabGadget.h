#pragma once

#include "gadgetronmatlab_export.h"
#include "Gadget.h"
#include "Gadgetron.h"
#include "hoNDArray.h"
#include "ismrmrd.h"

#include "engine.h"     // Matlab Engine header

#include <stdio.h>
#include <stdlib.h>

#include <complex>

template <class T> class MatlabGadget :
    public Gadget2<T, hoNDArray< std::complex<float> > >
{
public:
    MatlabGadget(): Gadget2<T, hoNDArray< std::complex<float> > >()
    {
        // Open the Matlab Engine on the current host
        if (!(engine_ = engOpen(""))) {
            // TODO: error checking!
            GADGET_DEBUG1("Can't start MATLAB engine\n");
        }
    }

    ~MatlabGadget()
    {
        // Close the Matlab engine
        engClose(engine_);
    }

protected:

    int process_config(ACE_Message_Block* mb)
    {

        path_        = this->get_string_value("path");
        classname_   = this->get_string_value("matlab_classname");

        GADGET_DEBUG2("MATLAB Path    : %s\n", path_.get()->c_str());
        GADGET_DEBUG2("MATLAB Class Name : %s\n", classname_.get()->c_str());

        // Set up buffer for catching Matlab output
#define BUFSIZE 4096
        char buffer[BUFSIZE] = "\0";
        engOutputBuffer(engine_, buffer, BUFSIZE);

        // Add +ismrmrd package to Matlab's path
        // TODO: this should be in user's Matlab path NOT HERE
        engEvalString(engine_, "addpath(strcat(getenv('ISMRMRD_HOME'), '/matlab'));");
        engEvalString(engine_, "addpath(strcat(getenv('GADGETRON_HOME'), '/matlab'));");

        // Check that we found ismrmrd package
        engEvalString(engine_, "which ismrmrd.AcquisitionHeader");

        // add user specified path for this gadget
        std::string add_user_path("addpath(");
        add_user_path += *path_.get() + ");";
        printf("%s\n", add_user_path.c_str());
        engEvalString(engine_, add_user_path.c_str());

        // Check that we found the class
        std::string which_classname("which ");
        which_classname += *classname_.get() + ";";
        engEvalString(engine_, which_classname.c_str());

        // Instantiate the Matlab gadget object from the user specified class
        std::string instantiate_gadget("matgadget = ");
        instantiate_gadget += *classname_.get() + "();";
        engEvalString(engine_, instantiate_gadget.c_str());

        // Call matlabgadget.config with the XML header
        std::string xmlConfig = std::string(mb->rd_ptr());
        mxArray *xmlhdr = mxCreateString(xmlConfig.c_str());
        engPutVariable(engine_, "xmlhdr", xmlhdr);
        engEvalString(engine_, "matgadget.config(xmlhdr);");

        GADGET_DEBUG2("%s", buffer);

	mxDestroyArray(xmlhdr);

        return GADGET_OK;
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
