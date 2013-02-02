function accumulator(head, data)
  global accumulator_calling_gadget;
  global accumulator_next_gadget;
  global accumulator_buffer;
  global accumulator_center_line;

   line_offset = bitshift(size(accumulator_buffer,2),-1) - accumulator_center_line;
   %size(accumulator_buffer(:,head.idx.kspace_encode_step_1+line_offset+1,1,:))
   %size(reshape(data,[size(data,1),1,1,size(data,2))])
   %size(data)
   
   accumulator_buffer(:,head.idx.kspace_encode_step_1+line_offset+1,1,:) = reshape(data,size(data,1),1,1,size(data,2));
   
  if (bitand(head.flags, bitshift(1,7)) > 0),


    img = ismrm_transform_kspace_to_image(accumulator_buffer,[1,2,3]);
    img = sqrt(sum(abs(img).^2,4));

    img_head = struct();
  
    img_head.version = h.version;
    img_head.flags = 0;
    img_head.measurement_uid = h.measurement_uid;
    img_head.matrix_size = [size(img,1),size(img,2),size(img,3)];
    img_head.channels = 1; 
    GadgetronReturnIsmrmrdImage(accumulator_next_gadget, head, single(img));
  end
end
