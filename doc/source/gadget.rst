Writing a Gadget
================

A Gadget is a :cpp:class:`Node<Gadgetron::Core::Node>` in the Gadgetron chain, 
which processes data comming in through an :cpp:class:`GenericInputChannel<Gadgetron::Core::GenericInputChannel>` and
sends the processed data to the next :cpp:class:`Node<Gadgetron::Core::Node>` in the chain using an
 :cpp:class:`OutputChannel<Gadgetron::Core::OutputChannel>`.

The simplest Gadgets to write are :cpp:class:`PureGadget<Gadgetron::Core::PureGadget>` and
:cpp:class:`ChannelGadget<Gadgetron::Core::ChannelGadget>`.

PureGadget
----------

A :cpp:class:`PureGadget<Gadgetron::Core::PureGadget>` is a Gadget which processes Messages one at a time,
and holds no state. Examples could be a Gadget which removes oversampling on :cpp:class:`Acquisitions<Gadgetron::Core::Acquisition>`,
or one which takes an :cpp:class:`Image<Gadgetron::Core::Image>` and performs autoscaling.

A PureGadget inhertis from :cpp:class:`PureGadget\<OUTPUT,INPUT\><Gadgetron::Core::PureGadget>`,
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

Channels are ranges, meaning they can be used directly with for loops and with algorithm from standard library, such as
:cpp:func:`std::transform` and :cpp:func:`std::accumulate`.

.. code-block: cpp

    void process(InputChannel<Acquisition>& in, OutputChannel& out) override {
        for (auto acquisition : in ) {

            auto& header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
            auto& data = std::get<hoNDArray<std::complex<float>>>(acquisition);
            auto& trajectory = std::get<optional<hoNDArray<float>>>(acquisition);
            //Gather acquisitions here
        }
    }

Or if you're using C++17, this would be

.. code-block: cpp

    void process(InputChannel<Acquisition>& in, OutputChannel& out) override {
        for (auto [header, data, trajectory] : in ) {
            //Gather acquisitions here
        }
    }





