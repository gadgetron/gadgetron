Tutorial: Writing a simple reconstruction Gadget
================================================

We wish to create a simple Gadget, which takes a buffer of fully sampled kspace data, and reconstructs them into a complex image.

The simplest form of Gadget is called a PureGadget, and creates one output for every input.

In order to create a PureGadget, we start with the following structure:

.. code-block: cpp

