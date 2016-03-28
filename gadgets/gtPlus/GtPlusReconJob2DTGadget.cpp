
#include "GtPlusReconJob2DTGadget.h"
#include "GtPlusGadgetOpenMP.h"

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

GtPlusReconJob2DTGadget::GtPlusReconJob2DTGadget()
{
    debugFolder_ = "DebugOutput";

    performTiming_ = true;

    verboseMode_ = false;

    gt_timer1_.set_timing_in_destruction(false);
    gt_timer2_.set_timing_in_destruction(false);
    gt_timer3_.set_timing_in_destruction(false);

    process_config_called_ = false;

    Gadgetron::prepOpenMP();
}

GtPlusReconJob2DTGadget::~GtPlusReconJob2DTGadget()
{

}

bool GtPlusReconJob2DTGadget::readParameters()
{
    try
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "------> GtPlusReconJob2DTGadget parameters <------");

        debugFolder_ = debugFolder.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "debugFolder_ is " << debugFolder_);

        performTiming_ = performTiming.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "performTiming_ is " << performTiming_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusReconJob2DTGadget::readParameters() ... ");
        return false;
    }

    return true;
}

int GtPlusReconJob2DTGadget::process_config(ACE_Message_Block* mb)
{
    // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
    //   0  1  2   3    4   5    6     7  8   9

    verboseMode_ = verboseMode.value();

    // read in parameters from the xml
    GADGET_CHECK_RETURN(this->readParameters(), GADGET_FAIL);

    // generate the destination folder
    if ( !debugFolder_.empty() )
    {
        Gadgetron::getDebugFolderPath(debugFolder_, debugFolder_fullPath_, verboseMode_);
    }
    else
    {
        GDEBUG_STREAM("GtPlusRecon, debugFolder is not set ...");
    }

    return GADGET_OK;
}

int GtPlusReconJob2DTGadget::process(Gadgetron::GadgetContainerMessage< int >* m1, Gadgetron::GadgetContainerMessage< GtPlusReconJobTypeCPFL > * m2)
{
    // because the parameter configuration will not be sent, we need to call process_config explicitly
    if ( !process_config_called_ )
    {
        GADGET_CHECK_RETURN( (this->process_config(m1)==0), GADGET_FAIL);
        process_config_called_ = true;
    }
    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusReconJob2DTGadget::process(...) starts ... ");

    int* jobID = m1->getObjectPtr();
    GDEBUG_CONDITION_STREAM(verboseMode_, "--> arriving job : " << *jobID << " ... ");

    GtPlusReconJobTypeCPFL* job = m2->getObjectPtr();
    GDEBUG_CONDITION_STREAM(verboseMode_, "    job array size : [ " << job->kspace.get_size(0) << " " 
                                                                 << job->kspace.get_size(1) << " " 
                                                                 << job->kspace.get_size(2) << " " 
                                                                 << job->kspace.get_size(3) << " ] ... ");

    // set the worker
    worker_grappa_.performTiming_ = performTiming_;
    if ( !debugFolder_fullPath_.empty() ) worker_grappa_.debugFolder_ = debugFolder_fullPath_;

    worker_noacceleration_.performTiming_ = performTiming_;
    if ( !debugFolder_fullPath_.empty() ) worker_noacceleration_.debugFolder_ = debugFolder_fullPath_;

    worker_spirit_.performTiming_ = performTiming_;
    if ( !debugFolder_fullPath_.empty() ) worker_spirit_.debugFolder_ = debugFolder_fullPath_;

    if ( verboseMode_ )
    {
        job->workOrder2DT.print(std::cout);
    }

    bool succeed = true;
    if ( performTiming_ ) { gt_timer1_.start("Recon 2DT job ... "); }

    succeed = worker_spirit_.performUnwarppingImpl(*job);

    if ( performTiming_ ) { gt_timer1_.stop(); }

    // export the results
    if ( !debugFolder_fullPath_.empty() )
    {
        std::ostringstream ostr;
        ostr << "ReconJob2DT_ID" << *jobID;

        hoNDArray< std::complex<float> > res = job->res;
        res.squeeze();
        if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_fullPath_+ostr.str()); }
    }

    // clean the kspace and ker and coil map
    job->kspace.clear();
    job->ker.clear();
    if ( !job->workOrder2DT.coilMap_ ) job->workOrder2DT.coilMap_->clear();

    if ( !succeed )
    {
        job->complexIm.clear();
        job->res.clear();
    }

    // send out the results
    GADGET_CHECK_RETURN(this->sendOutJob(*jobID, job), GADGET_FAIL);

    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusReconJob2DTGadget::process(...) ends ... ");

    m1->release();

    return GADGET_OK;
}

bool GtPlusReconJob2DTGadget::
sendOutJob(int jobID, GtPlusReconJobTypeCPFL* job)
{
    try
    {
      GDEBUG("GtPlusReconJob2DTGadget sendOutJob ...\n");

        if (!this->controller_)
        {
	  GERROR("Cannot return result to controller, no controller set\n");
	  return false;
        }

        GadgetContainerMessage<GadgetMessageIdentifier>* mb =
            new GadgetContainerMessage<GadgetMessageIdentifier>();

        mb->getObjectPtr()->id = GADGET_MESSAGE_CLOUD_JOB;

        GadgetContainerMessage<int>* m1 = new GadgetContainerMessage<int>();
        *(m1->getObjectPtr()) = jobID;

        GadgetContainerMessage<GtPlusReconJobTypeCPFL>* m2 = new GadgetContainerMessage<GtPlusReconJobTypeCPFL>();

        *(m2->getObjectPtr()) = *job;

        m1->cont(m2);
        mb->cont(m1);

        int ret =  this->controller_->output_ready(mb);
        if (ret < 0)
        {
            GDEBUG("Failed to return GtPlusReconJob2DTGadget job massage to controller\n");
            return false;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusReconJob2DTGadget::sendOutJob(...) ... ");
        return false;
    }

    return true;
}

GADGET_FACTORY_DECLARE(GtPlusReconJob2DTGadget)

}
