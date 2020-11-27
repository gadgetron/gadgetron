/** \file   GenericReconCartesianGrappaAIGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian grappa ai reconstruction, working on the IsmrmrdReconData.
    \author Hui Xue
*/

#pragma once

#include "GenericReconCartesianGrappaGadget.h"
#include "python_toolbox.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconCartesianGrappaAIGadget : public GenericReconCartesianGrappaGadget
    {
    public:
        GADGET_DECLARE(GenericReconCartesianGrappaAIGadget);

        typedef GenericReconCartesianGrappaGadget BaseClass;
        typedef typename BaseClass::ReconObjType ReconObjType;
        typedef std::complex<float> T;

        GenericReconCartesianGrappaAIGadget();
        ~GenericReconCartesianGrappaAIGadget();

    protected:

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        /// [Nor1 Sor1 SLC]
        std::vector< std::vector< hoNDArray<T> > > kernels_;
        std::vector< std::vector<boost::python::object> > models_;
        /// [RO E1 E2 1 N S SLC]
        std::vector<IsmrmrdImageArray> recon_res_grappa_ai_;

        // gadgetron home
        std::string gt_home_;

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);

        // calibration, if only one dst channel is prescribed, the GrappaOne is used
        virtual void perform_calib(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // unwrapping or coil combination
        virtual void perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        virtual int close(unsigned long flags);
    };
}
