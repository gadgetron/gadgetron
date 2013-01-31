function my_recon_function(head, data)
  global my_calling_gadget;
  global my_next_gadget;
  f = hamming(size(data,1));
  f = repmat(f,1,size(data,2));
  data = data .* f;
  GadgetronReturnIsmrmrdAcquisition(my_next_gadget, head, data);
end
