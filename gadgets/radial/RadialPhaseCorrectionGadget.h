#pragma once

#include "Gadget.h"
#include "hoNDArray.h"

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <complex>
#include <queue>
#include <map>

namespace Gadgetron {

  class RadialPhaseCorrectionGadget :
    public Gadget1<mrd::Acquisition>
  {
  public:
    RadialPhaseCorrectionGadget();
    ~RadialPhaseCorrectionGadget() {};

  protected:
    GADGET_PROPERTY_LIMITS(mode,int, "Radial mode", 3, GadgetPropertyLimitsEnumeration, 2,3);
    GADGET_PROPERTY(order,int,"Order of polynomial fit", 6);
    GADGET_PROPERTY(profiles, int, "Number of profiles to estimate fit", 500);

    virtual int process_config(const mrd::Header& header);

    virtual int process( GadgetContainerMessage<mrd::Acquisition> *m1);

    unsigned int mode_;
    unsigned int order_;
    unsigned int profiles_;
    unsigned int profiles_counter_;
    int slices_;
    int sets_;
    int channels_;

    using AcquisitionMessage = GadgetContainerMessage<mrd::Acquisition>;

    std::vector<bool> fit_calculated;
    std::vector<double> polyfit;
    std::map<unsigned int, std::queue<std::unique_ptr<AcquisitionMessage>>> profiles_queue;

//    boost::shared_array<bool> fit_calculated_;
//    boost::shared_array<double> polyfit_;
//    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > profiles_queue_;

  private:
    double get_projection_angle( unsigned int profile_idx );
    void phase_correct( GadgetContainerMessage<mrd::Acquisition>* );
  };
}
