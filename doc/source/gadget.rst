Writing a Gadget
================

A Gadget is a :cpp:class:`Node<Gadgetron::Core::Node>` in the Gadgetron chain, 
which processes data coming in through an :cpp:class:`GenericInputChannel<Gadgetron::Core::GenericInputChannel>` and
sends the processed data to the next :cpp:class:`Node<Gadgetron::Core::Node>` in the chain using an
 :cpp:class:`OutputChannel<Gadgetron::Core::OutputChannel>`.

The simplest Gadgets to write are :cpp:class:`PureGadget<Gadgetron::Core::PureGadget>` and
:cpp:class:`ChannelGadget<Gadgetron::Core::ChannelGadget>`.

PureGadget
----------

A :cpp:class:`PureGadget<Gadgetron::Core::PureGadget>` is a Gadget which processes Messages one at a time,
and holds no state. Examples could be a Gadget which removes oversampling on :cpp:class:`Acquisitions<Gadgetron::Core::Acquisition>`,
or one which takes an :cpp:class:`Image<Gadgetron::Core::Image>` and performs autoscaling.

A PureGadget inheritss from :cpp:class:`PureGadget\<OUTPUT,INPUT\><Gadgetron::Core::PureGadget>`,
where OUTPUT and INPUT are the output type and input type of the Gadget.

**AutoScaleGadget.h**

.. literalinclude:: ../../gadgets/mri_core/AutoScaleGadget.h
    :language: cpp

The **NODE_PROPERTY** macro defines a variable on the AutoScaleGadget which can be set from the XML file defining the chain.

**AutoScaleGadget.cpp**

.. literalinclude:: ../../gadgets/mri_core/AutoScaleGadget.cpp
    :language: cpp

Note the **GADGETRON_GADGET_EXPORT** declaration, which produces the code causing the AutoScaleGadget to be loadable by Gadgetron.

ChannelGadget
-------------

:cpp:class:`PureGadget<Gadgetron::Core::PureGadget>` can't hold any state between different messages and must send
one message per input.
This makes it easier to reason about and implement, but is also limiting in cases where we want to accumulate multiple
messages for processing. In this case, :cpp:class:`ChannelGadget<Gadgetron::Core::ChannelGadget>` should be used.
If we want to create a Gadget which takes several :cpp:class:`Acquisitions<Gadgetron::Core::Acquisition>` and
reconstruct them, we inherit from :cpp:class:`ChannelGadget\<Acquisition\><Gadgetron::Core::ChannelGadget>`.

.. code-block:: cpp

    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    class SimpleRecon : public ChannelGadget<Acquisition>{
        public:
            void process(InputChannel<Acquisition>& in, OutputChannel& out) override {
                ...
            }
    }

We can take the messages from the channel either by calling :cpp:func:`InputChannel::pop()` directly, or by using it in
 a for loop.
Channels are ranges, meaning they can be used directly with for loops and with algorithm from standard library, such as
:cpp:func:`std::transform` and :cpp:func:`std::accumulate`.

.. code-block:: cpp

    void process(InputChannel<Acquisition>& in, OutputChannel& out) override {
        for (auto acquisition : in ) {

            auto& header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
            auto& data = std::get<hoNDArray<std::complex<float>>>(acquisition);
            auto& trajectory = std::get<optional<hoNDArray<float>>>(acquisition);
            //Gather acquisitions here
        }
    }

Or if you're using C++17, this would be

.. code-block:: cpp

    void process(InputChannel<Acquisition>& in, OutputChannel& out) override {
        for (auto [header, data, trajectory] : in ) {
            //Gather acquisitions here
        }
    }

We want to gather acquisitions until we have enough for a (possibly undersampled) image. The AcquisitionHeader has the
 ISMRMRD::_ACQ_LAST_IN_ENCODE_STEP1 flag which we can use as a trigger. By importing channel_algorithms.h, we can write

.. code-block:: cpp

    void process(InputChannel<Acquisition>& in, OutputChannel& out) override {

        auto split_condition = [](auto& message){
          return std::get<ISMRMRD::AcquisitionHeader>(message).isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1);
        };

        for (auto acquisitions : buffer(in,split_condition)) {
            for (auto [header, data, trajectory] : acquisitions ) {
            //Gather acquisitions here
            }
    }


.. code-block:: cpp

    #include <gadgetron/Gadget.h>
    #include <gadgetron/hoNDFFT.h>
    #include <gadgetron/mri_core_utility.h>
    #include <gadgetron/ChannelAlgorithms.h>
    #include <gadgetron/log.h>
    #include <gadgetron/mri_core_coil_map_estimation.h>
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Algorithm;

    class SimpleRecon : public ChannelGadget<Acquisition> {

        public:
            SimpleRecon(const Context& context, const GadgetProperties& params) : ChannelGadget<Acquisition>(params), header{context.header} {

            }

            void process(InputChannel<Acquisition>& in, OutputChannel& out){

                auto recon_size = header.encoding[0].encodedSpace.matrixSize;

                ISMRMRD::AcquisitionHeader saved_header;

                auto split_condition = [](auto& message){
                return std::get<ISMRMRD::AcquisitionHeader>(message).isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1);
                };

                for (auto acquisitions : buffer(in,split_condition)) {

                   auto data = hoNDArray<std::complex<float>>(recon_size.x,recon_size.y,recon_size.z,header.acquisitionSystemInformation->receiverChannels.get());
                   for ( auto [acq_header, acq_data, trajectories] : acquisitions){
                        saved_header = acq_header;
                        data(slice,acq_header.idx.kspace_encode_step_1,0,slice) = acq_data;
                    }

                    hoNDFFT<float>::instance()->fft2c(data);

                    auto coil_map = coil_map_Inati(data);
                    data = coil_combine(data,coil_map,3);

                    auto image_header = image_header_from_acquisition(saved_header,header,data);

                    out.push(image_header,data);
                }
            }
        private:
            const ISMRMRD::IsmrmrdHeader header;
    };

    GADGETRON_GADGET_EXPORT(SimpleRecon)

   