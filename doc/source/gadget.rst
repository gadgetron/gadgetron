Writing a Gadget
================

A Gadget is a :cpp:class:`Node<Gadgetron::Core::Node>` in the Gadgetron chain, 
which processes data comming in through an :cpp:class:`InputChannel<Gadgetron::Core::InputChannel>` and 
sends the processed data to the next :cpp:class:`Node<Gadgetron::Core::Node>` in the chain using an :cpp:class:`OutputChannel<Gadgetron::Core::OutputChannel>`.

The simplest Gadgets to write are :cpp:class:`(Typed)PureGadget<Gadgetron::Core::TypedPureGadget>` and 
:cpp:class:`(Typed)ChannelGadget<Gadgetron::Core::TypedChannelGadget>`. 

PureGadget
----------

A :cpp:class:`PureGadget<Gadgetron::Core::TypedPureGadget>` is a Gadget which processes Messages one at a time,
and holds no state. Examples could be a Gadget which removes oversampling on :cpp:class:`Acquisitions<Gadgetron::Core::Acquisition>`,
or one which takes an :cpp:class:`Image<Gadgetron::Core::Image>` and performs autoscaling.

A TypedPureGadget inhertis from :cpp:class:`TypedPureGadget\<OUTPUT,INPUT\><Gadgetron::Core::TypedPureGadget>`,
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



