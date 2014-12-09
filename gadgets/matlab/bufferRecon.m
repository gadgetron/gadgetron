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
            nx = g.xml.encoding.encodedSpace.matrixSize.x;
            ny = g.xml.encoding.encodedSpace.matrixSize.y;
            nz = g.xml.encoding.encodedSpace.matrixSize.z;
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
            
            data = recon_data{1}.data; %Assume only one encoding space
            head = recon_data{1}.headers{1}; %Just get header from first trajectory
            
            
            % At the end of the acquisition, reconstruct the slice
            img_head = ismrmrd.ImageHeader;

            %img_head.slice = head.idx.slice;
            % set the matrix size
            % set one element at a time to not break the type (uint16) of matrix_size
            img_head.matrix_size(1) = g.img_size(1); % nx
            img_head.matrix_size(2) = g.img_size(2); % ny
            img_head.matrix_size(3) = g.img_size(3); % nz
            
            img_head.position = head.position;
            img_head.read_dir = head.read_dir;
            img_head.phase_dir = head.phase_dir;
            img_head.slice_dir = head.slice_dir;
            img_head.patient_table_position = head.patient_table_position;
            img_head.acquisition_time_stamp = head.acquisition_time_stamp;
            img_head.image_index = g.image_num;
            img_head.image_series_index = g.series_num;
            
            img_data = squeeze(data);
            img_data = fftshift(ifftn(fftshift(img_data)));
            imagesc(abs(img_data(:,:,1,1))); axis image; axis square;
            pause(2)
            close()
            
            disp(size(img_data));
            g.putImageQ(img_head, img_data);
            %fprintf('Put on Queue %d, type = %d\n',length(g.Q),g.Q{1}.type);
            
            
            
        end
        
    end
end
