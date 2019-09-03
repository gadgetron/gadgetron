Type matching
=============

:cpp:class:`PureGadget<Gadgetron::Core::PureGadget>`, :cpp:class:`ChannelGadget<Gadgetron::Core::ChannelGadget>` as well as several other classes in Gadgetron have a template arguement defining what type of data they accept.
In the simplest case, this is simply a list of types. For instance 

.. code-block:: cpp

    ChannelGadget<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float>>>

will match any message starting with a :cpp:class:`AcquisitionHeader` followed by a :cpp:class:`hoNDArray<std::complex<float>>`.
If the message has more parts, these will simply be discarded. Also note that this is equivalent to 

.. code-block:: cpp

    ChannelGadget<tuple<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float>>>>

Optional
--------

If we want to include an element that will only appear some times, we can define it as optional. For instance, acquisitions can have a trajectory attached to them. This would look like 

.. code-block:: cpp

    ChannelGadget<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float>>, optional<hoNDArray<float>>>

and in fact Types.h defines :cpp:class:`Acquisition<Gadgetron::Core::Acquisition>` as 

.. code-block:: cpp

    using Acquisition = tuple<ISMRMRD::AcquisitionHeader,  hoNDArray<std::complex<float>>,optional<hoNDArray<float>>>;




Variant
-------

What if you need to create a ChannelGadget that accepts multiple types? For instance, one which receives both Acquisition and Waveform. In this case we can use a variant.

.. code_block:: cpp

   ChannelGadget<variant<Acquisition,Waveform>> 

In order to work with the data, you call :cpp:func:`Core::visit<Gadgetron::Core::visit>`, which is modelled from `std::visit <https://en.cppreference.com/w/cpp/utility/variant/visit>`_.

For instance, a toy example which counts the number of data points in all waveforms and acquisitions could look like

.. code_block:: cpp

   void process(ChannelGadget<Acquisition,Waveform>& in, OutputChannel& out>){
       
       for (auto acquisition_or_waveform : in){
           
           size_t counts = 0;
           Core::visit( [counts&](auto& val){
               auto& data = std::get<1>(val); //Data the second arguement for both acquisitons and waveforms
               counts += data.size();
           },
           acquisition_or_waveform);
       }
   }





