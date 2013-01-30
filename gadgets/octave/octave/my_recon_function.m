function my_recon_function(head, data)
  global my_calling_gadget;
  global my_next_gadget;

  %fprintf("My recon function called with line = %d\n", head.idx.kspace_encode_step_1);
  %fprintf("My recon function called with line = %d\n", head.scan_counter);

  %GadgetronReturnIsmrmrd(my_calling_gadget);
  %plot(abs(data(:,1)));

  f = hamming(size(data,1));
  f = repmat(f,1,size(data,2));
  sum(abs(data(:)).^2)
  data = data .* f;
  sum(abs(data(:)).^2)
  GadgetronReturnIsmrmrdAcquisition(my_next_gadget, head, data);
end
