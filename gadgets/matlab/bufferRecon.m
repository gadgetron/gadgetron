classdef bufferRecon < handle & BaseBufferGadget
    
    properties
        
        image_num;
        series_num;
        center_line;
        img_size;
        
    end
    
    methods
        
        function g = config(g)
            
            fprintf('The resonance frequency is %d\n', g.xml.experimentalConditions.H1resonanceFrequency_Hz);
            nx = g.xml.encoding.reconSpace.matrixSize.x;
            ny = g.xml.encoding.reconSpace.matrixSize.y;
            nz = g.xml.encoding.reconSpace.matrixSize.z;
            % for 2D sequences the number of getZ breaks
            %try
            %catch
            
            % nz =1
            %end
            % the number of receiver channels is optional
            try
                % this is the only cast from java.lang.Integer that works in Matlab
                nc = g.xml.acquisitionSystemInformation.receiverChannels;
            catch
                nc = 1;
            end
            % the number of slices is optional
            try
                ns = g.xml.encoding.encodingLimits.slice.maximum + 1;
            catch
                ns = 1;
            end
            
            g.center_line = g.xml.encoding.encodingLimits.kspace_encoding_step_1.center;
            g.img_size = [nx ny nz];
            g.image_num = 0;   % todo this needs to be static or global...
            g.series_num = 0;  % todo this needs to be static or global...
        end
        
        function g = process(g, recon_data)
            disp('Processing')
            % stuff the line
            recon_data.data
            data = recon_data(1).data; %Assume only one encoding space
            head = data.headers; %Just get header from first trajectory
            
            img_data = squeeze(data.data);
            img_data = fftshift(fftshift(ifft2(ifftshift(ifftshift(img_data,1),2)),1),2);
            img_data = squeeze(sqrt(sum(img_data.^2,3)));
            %img_data = img_data./sqrt(size(img_data,1)*size(img_data,2));
            
            % At the end of the acquisition, reconstruct the slice
            img_head = ismrmrd.ImageHeader;


            % set the matrix size
            % set one element at a time to not break the type (uint16) of matrix_size
            img_head.matrix_size(1) = size(img_data,1);
            img_head.matrix_size(2) = size(img_data,2);
            img_head.matrix_size(3) = size(img_data,3);

            img_head.position = head.position(:,1);
            img_head.read_dir = head.read_dir(:,1);
            img_head.phase_dir = head.phase_dir(:,1);
            img_head.slice_dir = head.slice_dir(:,1);
            img_head.patient_table_position = head.patient_table_position(:,1);
            img_head.acquisition_time_stamp = head.acquisition_time_stamp(1);
            img_head.image_index = g.image_num;
            img_head.image_series_index = g.series_num;
            img_head.channels = 1;
            
           
            %imagesc(abs(img_data(:,:,1,1))); axis image; axis square;
            %pause(2)
            %close()
            
            %disp(size(img_data));
            
            g.putImageQ(img_head, img_data);
            %fprintf('Put on Queue %d, type = %d\n',length(g.Q),g.Q{1}.type);
            
            
            
        end
        
    end
end
