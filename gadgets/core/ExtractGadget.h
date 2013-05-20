#ifndef EXTRACTGADGET_H_
#define EXTRACTGADGET_H_

#include "Gadget.h"
#include "gadgetron_core_export.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"

#include <complex>

#define MAX_UNSIGNED_SHORT_IMAGE_VALUE

//Extract flags
#define GADGET_EXTRACT_NONE                   (0)      //0
#define GADGET_EXTRACT_MAGNITUDE              (1 << 0) //1
#define GADGET_EXTRACT_REAL                   (1 << 1) //2
#define GADGET_EXTRACT_IMAG                   (1 << 2) //4
#define GADGET_EXTRACT_PHASE                  (1 << 3) //8
#define GADGET_EXTRACT_MAX                    (1 << 4) //16

namespace Gadgetron{

class EXPORTGADGETSCORE ExtractGadget:
public Gadget2<ISMRMRD::ImageHeader,hoNDArray< std::complex<float> > >
{

public:
  GADGET_DECLARE(ExtractGadget);

  ExtractGadget();
  virtual ~ExtractGadget();

  void set_extract_mask(unsigned short mask) {
	  extract_mask_ = mask;
  }

  bool extract_magnitude() {
	  return (extract_mask_ & GADGET_EXTRACT_MAGNITUDE);
  }

  bool extract_real() {
	  return (extract_mask_ & GADGET_EXTRACT_REAL);
  }

  bool extract_imag() {
	  return (extract_mask_ & GADGET_EXTRACT_IMAG);
  }

  bool extract_phase() {
	  return (extract_mask_ & GADGET_EXTRACT_PHASE);
  }

protected:
	virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

	unsigned short extract_mask_;
};
}

#endif /* EXTRACTGADGET_H_ */
