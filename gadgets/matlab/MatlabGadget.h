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
        reffunc_     = this->get_string_value("gadget_reference_function");
        datafunc_    = this->get_string_value("input_function");
        configfunc_  = this->get_string_value("config_function");

        GADGET_DEBUG2("MATLAB Ref Function    : %s\n", reffunc_.get()->c_str());
        GADGET_DEBUG2("MATLAB Data Function   : %s\n", datafunc_.get()->c_str());
        GADGET_DEBUG2("MATLAB Config Function : %s\n", configfunc_.get()->c_str());

        // Parse ISMRMRD XML header
        //boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(string(mb->rd_ptr()));
        std::string xmlConfig = std::string(mb->rd_ptr());
        mxArray *xmlhdr = mxCreateString(xmlConfig.c_str());
        engPutVariable(engine_, "xmlhdr", xmlhdr);

        char buffer[1025] = "\0";
        engOutputBuffer(engine_, buffer, 1024);

        engEvalString(engine_, "fprintf(xmlhdr)");

        printf("%s", buffer);

	mxDestroyArray(xmlhdr);

        return GADGET_OK;
    }

    boost::shared_ptr<std::string> path_;
    boost::shared_ptr<std::string> reffunc_;
    boost::shared_ptr<std::string> datafunc_;
    boost::shared_ptr<std::string> configfunc_;

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
