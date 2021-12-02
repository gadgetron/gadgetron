#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_radial_export.h"

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <queue>
#include <map>

namespace Gadgetron {

  class EXPORTGADGETS_RADIAL RadialPhaseCorrectionGadget :
    public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
  {
  public:
    GADGET_DECLARE(RadialPhaseCorrectionGadget);
    RadialPhaseCorrectionGadget();
    ~RadialPhaseCorrectionGadget() {};
    
  protected:
    GADGET_PROPERTY_LIMITS(mode,int, "Radial mode", 3, GadgetPropertyLimitsEnumeration, 2,3);
    GADGET_PROPERTY(order,int,"Order of polynomial fit", 6);
    GADGET_PROPERTY(profiles, int, "Number of profiles to estimate fit", 500);
    
    virtual int process_config( ACE_Message_Block *mb );
    
    virtual int process( GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
                         GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2);
    
    unsigned int mode_;
    unsigned int order_;
    unsigned int profiles_;
    unsigned int profiles_counter_;
    int slices_;
    int sets_;
    int channels_;

    using AcquisitionMessagePtr = GadgetContainerMessage<ISMRMRD::AcquisitionHeader>*;

    std::vector<bool> fit_calculated;
    std::vector<double> polyfit;
    std::map<unsigned int, std::queue<AcquisitionMessagePtr>> profiles_queue;

//    boost::shared_array<bool> fit_calculated_;
//    boost::shared_array<double> polyfit_;
//    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > profiles_queue_;

  private:
    double get_projection_angle( unsigned int profile_idx );
    void phase_correct( GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* );
  };
}
