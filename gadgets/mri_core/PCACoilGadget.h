#ifndef PCACOILGADGET_H_
#define PCACOILGADGET_H_

#include "gadgetron_mricore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "hoNDKLT.h"
#include "ismrmrd/ismrmrd.h"

#include <complex>
#include <map>

namespace Gadgetron {

  class EXPORTGADGETSMRICORE PCACoilGadget :
    public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
  {
    typedef Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > > inherited;
  public:
    GADGET_DECLARE(PCACoilGadget);

    PCACoilGadget();
    virtual ~PCACoilGadget();

  protected:
    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

  private:
    GADGET_PROPERTY(uncombined_channels_by_name, std::string, "List of comma separated channels by name", "");
    GADGET_PROPERTY(present_uncombined_channels, int, "Number of uncombined channels found", 0);

    std::vector<unsigned int> uncombined_channels_;
    
    //Map containing buffers, one for each location
    std::map< int, std::vector< ACE_Message_Block* > > buffer_;

    //Keep track of whether we are buffering for a particular location
    std::map< int, bool> buffering_mode_;

    //Map for storing PCA coefficients for each location
    std::map<int, hoNDKLT<std::complex<float> >* > pca_coefficients_;

    int max_buffered_profiles_;
    int samples_to_use_;
  };
}

#endif /* PCACOILGADGET_H_ */
