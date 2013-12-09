/** \file   GtPlusRecon2DTGadgetCloud.h
    \brief  This is the gateway gadget for the dual layer GtPlus cloud.
            For every incoming k-space data package, it is sent to a first layer gadget.
            If a data package was not processed successfully and results were not returned to this gadget,
            the reconstruction will be performed locally.

            Ref to: 

            Hui Xue, Souheil Inati, Thomas Sangild Sorensen, Peter Kellman, Michael S. Hansen. 
            Distributed MRI Reconstruction using Gadgetron based Cloud Computing. Submitted to
            Magenetic Resonance in Medicine on Dec 2013.

    \author Hui Xue
*/

#pragma once

#include "GtPlusRecon2DTGadget.h"
#include "GadgetCloudController.h"
#include "GadgetCloudJobMessageReadWrite.h"
#include "GtPlusRecon2DTCloudPackage.h"

namespace Gadgetron
{

class EXPORTGTPLUSGADGET GtPlusRecon2DTGadgetCloud : public GtPlusRecon2DTGadget
{
public:
    GADGET_DECLARE(GtPlusRecon2DTGadgetCloud);

    typedef GtPlusRecon2DTGadget BaseClass;

    typedef BaseClass::ValueType ValueType;
    typedef BaseClass::WorkOrderType WorkOrderType;
    typedef BaseClass::WorkOrder2DTType WorkOrder2DTType;
    typedef BaseClass::DimensionRecordType DimensionRecordType;

    typedef GtPlusRecon2DTCloudPackage<ValueType> CloudPackageType;

    typedef Gadgetron::GadgetCloudController< CloudPackageType > GTCloudControllerType;

    GtPlusRecon2DTGadgetCloud();
    ~GtPlusRecon2DTGadgetCloud();

    virtual int close(unsigned long flags);

    std::vector<CloudPackageType> packages_sent_;
    std::vector<CloudPackageType> packages_received_;

    // indicate whether the results of all sent packages have been passed to next gadget or not
    std::vector< std::pair<unsigned int, bool> >  packages_passed_to_next_gadget_;

    // store the image headers for every incoming package
    std::vector<GtPlusGadgetImageArray> image_headers_;

protected:

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(Gadgetron::GadgetContainerMessage< GtPlusGadgetImageArray >* m1, Gadgetron::GadgetContainerMessage< WorkOrderType > * m2);

    virtual bool processJob(CloudPackageType& jobSent, CloudPackageType& jobReceived);

    GTCloudControllerType controller_;

    unsigned int curr_node_;

    unsigned int num_of_jobs_;

    std::vector<GadgetMessageReader*> readers_;
    std::vector<GadgetMessageWriter*> writers_;

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer_2DT_cloud_;
};

class GtPlusRecon2DTGadgetCloudSender : public GadgetCloudJobProcessHandler< GtPlusRecon2DTCloudPackage< std::complex<float> > >
{
public:

    typedef std::pair<Gadgetron::gtPlus::ISMRMRDDIM, unsigned long long> DimensionRecordType;

    GtPlusRecon2DTGadgetCloudSender();
    virtual ~GtPlusRecon2DTGadgetCloudSender();

    virtual bool processJob(int jobID, GtPlusRecon2DTCloudPackage< std::complex<float> >& ajob);

    // pointer to the gadget
    GtPlusRecon2DTGadgetCloud* gadget_;
};

}
