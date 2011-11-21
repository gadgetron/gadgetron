#include "MRIAcquisitionReader.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include <complex>
#include "GadgetMessageInterface.h"

ACE_Message_Block* MRIAcquisitionReader::read(ACE_SOCK_Stream* sock)
{
  
  GadgetContainerMessage<GadgetMessageAcquisition>* m1 = 
    new GadgetContainerMessage<GadgetMessageAcquisition>();
  
  GadgetContainerMessage<hoNDArray< std::complex<float> > >* m2 = 
    new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
  
  m1->cont(m2);
  
  ssize_t recv_cnt = 0;
  if ((recv_cnt = 
       sock->recv_n (m1->getObjectPtr(), 
		     sizeof(GadgetMessageAcquisition))) <= 0) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) Unable to read Acq header\n")));
    
    m1->release();
    
    return 0;
  }

  std::vector<unsigned int> adims;
  adims.push_back(m1->getObjectPtr()->samples);
  adims.push_back(m1->getObjectPtr()->channels);

  static int counter = 0;
  if (counter < 100) {
  	std::cout << "Receiving" << counter++ << ", " << m1->getObjectPtr()->idx.line << std::endl;
  }

  if (!m2->getObjectPtr()->create(&adims)) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) Allocate sample data\n")));

    m1->release();

    return 0;
  }
  
  if ((recv_cnt = 
       sock->recv_n
       (m2->getObjectPtr()->get_data_ptr(), 
	sizeof(std::complex<float>)*m1->getObjectPtr()->samples*m1->getObjectPtr()->channels)) <= 0) {

    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) Unable to read Acq data\n")));

    m1->release();

    return 0;
  }

  return m1;
}

GADGETRON_READER_FACTORY_DECLARE(MRIAcquisitionReader)
