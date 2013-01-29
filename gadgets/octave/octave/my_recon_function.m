function my_recon_function(tmp)
  global my_calling_gadget;
  fprintf("My recon function called with %s", tmp);
  GadgetronReturnIsmrmrd(my_calling_gadget);
end
