/*
 * PCACoilGadget.h
 *
 *  Created on: Dec 13, 2011
 *      Author: Michael S. Hansen
 */

#ifndef PCACOILGADGET_H_
#define PCACOILGADGET_H_

#include "gadgetronmricore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd.h"

#include <complex>
#include <map>
namespace Gadgetron{
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
	  //Map containing buffers, one for each location
      std::map< int, std::vector< ACE_Message_Block* > > buffer_;

      //Keep track of whether we are buffering for a particular location
      std::map< int, bool> buffering_mode_;

      //Map for storing PCA coefficients for each location
      std::map<int, hoNDArray<std::complex<float> >* > pca_coefficients_;

      int max_buffered_profiles_;
      int samples_to_use_;
};
}
#endif /* PCACOILGADGET_H_ */
