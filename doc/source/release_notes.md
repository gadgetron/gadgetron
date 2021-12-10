# Release Notes

## 4.0

+ New [interface](gadget) for Gadgets, based on [Channels](_channels), which makes it significantly easier to write new Gadgets. Old style Gadget interface still supported and are 100% compatible with the new interface.
+ Branching chains now supported
+ Vastly better support for distributed computing, which now allows for any Gadget which receives and sends standard Gadgetron datatypes to be distributed across many Gadgets.
+ Improved error handling. Gadgetron should now produce more meaningful errors and single Gadgets can no longer crash the Gadgetron instance.
+ Removed ACE. A compatibility layer has been added to provide a standin for ACE_Message_Block. Importing ACE headers in any shape or form is no longer supported.
+ Added support for NumPy-like slicing for hoNDArrays.

