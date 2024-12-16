Type matching
=============

:cpp:class:`PureGadget<Gadgetron::Core::PureGadget>`, :cpp:class:`ChannelGadget<Gadgetron::Core::ChannelGadget>` as well as several other classes in Gadgetron have a template argument defining what type of data they accept.
In the simplest case, this is simply a list of types. For instance

.. code-block:: cpp

    ChannelGadget<mrd::Acquisition>

will match any message containing an :cpp:class:`Acquisition`.

Optional
--------

If we want to include an element that will only appear some times, we can define it as optional. For instance, acquisitions can have a trajectory attached to them. This would look like

.. code-block:: cpp

    ChannelGadget<mrd::Acquisition, std::optional<hoNDArray<float>>>

Variant
-------

What if you need to create a ChannelGadget that accepts multiple types? For instance, one which receives both Acquisition and Waveform. In this case we can use a variant.

.. code-block:: cpp

   ChannelGadget<std::variant<mrd::Acquisition, mrd::Waveform>>

In order to work with the data, you call :cpp:func:`std::visit <https://en.cppreference.com/w/cpp/utility/variant/visit>`_.

For instance, a toy example which counts the number of data points in all waveforms and acquisitions could look like

.. code-block:: cpp

   void process(ChannelGadget<Acquisition, Waveform>& in, OutputChannel& out){

       for (auto acquisition_or_waveform : in){

           size_t counts = 0;
           std::visit( [counts&](auto& val){
               auto& data = std::get<1>(val); //Data the second argument for both acquisitons and waveforms
               counts += data.size();
           },
           acquisition_or_waveform);
       }
   }
