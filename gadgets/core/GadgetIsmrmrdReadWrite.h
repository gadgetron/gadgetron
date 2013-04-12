#ifndef GADGETISMRMRDREADWRITE_H
#define GADGETISMRMRDREADWRITE_H

#include "Gadgetron.h"
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "GadgetMessageInterface.h"
#include "hoNDArray.h"
#include "url_encode.h"
#include "gadgetroncore_export.h"

#include <ismrmrd.h>

#ifndef EXCLUDE_ISMRMRD_XSD
#include "ismrmrd.hxx"
#endif

#include <ace/SOCK_Stream.h>
#include <ace/Task.h>
#include <complex>

namespace Gadgetron{

    class GadgetIsmrmrdAcquisitionMessageWriter : public GadgetMessageWriter
    {

    public:
        virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb)
        {
            GadgetContainerMessage<ISMRMRD::Acquisition>* acqmb =
                dynamic_cast< GadgetContainerMessage<ISMRMRD::Acquisition>* >(mb);

            if (!acqmb) {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), GadgetAcquisitionMessageWriter, invalid acquisition message objects")) );
                return -1;
            }

            ssize_t send_cnt = 0;

            GadgetMessageIdentifier id;
            id.id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;

            if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
                ACE_DEBUG ((LM_ERROR,
                    ACE_TEXT ("(%P|%t) Unable to send acquisition message identifier\n")));

                return -1;
            }

            if ((send_cnt = sock->send_n (&acqmb->getObjectPtr()->getHead(), sizeof(ISMRMRD::AcquisitionHeader))) <= 0) {
                ACE_DEBUG ((LM_ERROR,
                    ACE_TEXT ("(%P|%t) Unable to send acquisition header\n")));

                return -1;
            }

            unsigned long trajectory_elements = acqmb->getObjectPtr()->getHead().trajectory_dimensions*acqmb->getObjectPtr()->getHead().number_of_samples;
            unsigned long data_elements = acqmb->getObjectPtr()->getHead().active_channels*acqmb->getObjectPtr()->getHead().number_of_samples;

            if (trajectory_elements) {
                if ((send_cnt = sock->send_n (&acqmb->getObjectPtr()->getTraj()[0], sizeof(float)*trajectory_elements)) <= 0) {
                    ACE_DEBUG ((LM_ERROR,
                        ACE_TEXT ("(%P|%t) Unable to send acquisition trajectory elements\n")));

                    return -1;
                }
            }

            if (data_elements) {
                if ((send_cnt = sock->send_n (&acqmb->getObjectPtr()->getData()[0], 2*sizeof(float)*data_elements)) <= 0) {
                    ACE_DEBUG ((LM_ERROR,
                        ACE_TEXT ("(%P|%t) Unable to send acquisition data elements\n")));

                    return -1;
                }
            }

            return 0;
        }
    };

    /**
    Default implementation of GadgetMessageReader for IsmrmrdAcquisition messages

    */
    class EXPORTGADGETSCORE GadgetIsmrmrdAcquisitionMessageReader : public GadgetMessageReader
    {

    public:
        GADGETRON_READER_DECLARE(GadgetIsmrmrdAcquisitionMessageReader);

        virtual ACE_Message_Block* read(ACE_SOCK_Stream* stream)
        {

            GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1 =
                new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();

            GadgetContainerMessage<hoNDArray< std::complex<float> > >* m2 =
                new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

            m1->cont(m2);

            ssize_t recv_count = 0;

            if ((recv_count = stream->recv_n(m1->getObjectPtr(), sizeof(ISMRMRD::AcquisitionHeader))) <= 0) {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetIsmrmrdAcquisitionMessageReader, failed to read ISMRMRDACQ Header\n")) );
                m1->release();
                return 0;
            }

            if (m1->getObjectPtr()->trajectory_dimensions) {
                GadgetContainerMessage<hoNDArray< float > >* m3 =
                    new GadgetContainerMessage< hoNDArray< float > >();

                m2->cont(m3);

                std::vector<unsigned int> tdims;
                tdims.push_back(m1->getObjectPtr()->trajectory_dimensions);
                tdims.push_back(m1->getObjectPtr()->number_of_samples);

                try { m3->getObjectPtr()->create(&tdims);}
                catch (runtime_error &err){
                    GADGET_DEBUG_EXCEPTION(err,"(%P|%t) Allocate trajectory data\n");
                    m1->release();

                    return 0;
                }

                if ((recv_count =
                    stream->recv_n
                    (m3->getObjectPtr()->get_data_ptr(),
                    sizeof(float)*tdims[0]*tdims[1])) <= 0) {

                        ACE_DEBUG ((LM_ERROR,
                            ACE_TEXT ("(%P|%t) Unable to read trajectory data\n")));

                        m1->release();

                        return 0;
                }

            }

            std::vector<unsigned int> adims;
            adims.push_back(m1->getObjectPtr()->number_of_samples);
            adims.push_back(m1->getObjectPtr()->active_channels);

            try{ m2->getObjectPtr()->create(&adims); }
            catch (runtime_error &err ){
                GADGET_DEBUG_EXCEPTION(err,"(%P|%t) Allocate sample data\n")
                    m1->release();

                return 0;
            }

            if ((recv_count =
                stream->recv_n
                (m2->getObjectPtr()->get_data_ptr(),
                sizeof(std::complex<float>)*adims[0]*adims[1])) <= 0) {

                    ACE_DEBUG ((LM_ERROR,
                        ACE_TEXT ("(%P|%t) Unable to read Acq data\n")));

                    m1->release();

                    return 0;
            }

            return m1;
        }

    };

#ifndef EXCLUDE_ISMRMRD_XSD
    inline boost::shared_ptr<ISMRMRD::ismrmrdHeader> parseIsmrmrdXMLHeader(std::string xml) {
        char * gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");
        ACE_TCHAR schema_file_name[4096];
        ACE_OS::sprintf(schema_file_name, "%s/schema/ismrmrd.xsd", gadgetron_home);

        std::string tmp(schema_file_name);
        tmp = url_encode(tmp);
        ACE_OS_String::strncpy(schema_file_name,tmp.c_str(), 4096);

        xml_schema::properties props;
        props.schema_location (
            "http://www.ismrm.org/ISMRMRD",
            std::string (schema_file_name));


        std::istringstream str_stream(xml, std::stringstream::in);

        boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg;

        ACE_TCHAR port_no[1024];
        try {
            cfg = boost::shared_ptr<ISMRMRD::ismrmrdHeader>(ISMRMRD::ismrmrdHeader_ (str_stream,0,props));
        }  catch (const xml_schema::exception& e) {
            GADGET_DEBUG2("Failed to parse XML Parameters: %s\n", e.what());
        }

        return cfg;
    }
#endif
}
#endif //GADGETISMRMRDREADWRITE_H
