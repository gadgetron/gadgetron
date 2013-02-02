function downsample_2x(head, data)
  global downsample_2x_calling_gadget;
  global downsample_2x_next_gadget;

  data = ismrm_transform_kspace_to_image(data,1);
  data = data([1:bitshift(size(data,1),-1)]+bitshift(size(data,1),-2),:);
  data = ismrm_transform_image_to_kspace(data,1);
  head.number_of_samples = size(data,1);

  GadgetronReturnIsmrmrdAcquisition(downsample_2x_next_gadget, head, single(data));
end
