classdef bufferRecon < handle & BaseBufferGadget
    
    properties
        
        image_num;
        series_num;
        
    end
    
    methods
        
        function g = config(g)
            
            %fprintf('The resonance frequency is %d\n', g.xml.experimentalConditions.H1resonanceFrequency_Hz);            
            g.image_num = 0;   % todo this needs to be static or global...
            g.series_num = 0;  % todo this needs to be static or global...
            
        end
        
        function g = process(g, recon_data)
            fprintf('Processing')
            
            %recon_data.data
            
            data = recon_data(1).data; %Assume only one encoding space
            head = data.headers; %Just get header from first trajectory
            
            img_data = data.data;
            % fft along x
            img_data = fftshift(ifft(fftshift(img_data,1),[],1),1);
            % fft along y
            img_data = fftshift(ifft(fftshift(img_data,2),[],2),2);
            if size(img_data,3) > 1
                % fft along y
                img_data = fftshift(ifft(fftshift(img_data,3),[],3),3);
            end
            
            % sqrt of sum of squares over channels
            img_data = sqrt(sum(abs(img_data).^2,4));
            
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
            
            %center = floor(size(img_data,3)/2)+1;
            %imagesc(abs(img_data(:,:,center,1,1,1,1))); axis image; axis square;
            %pause(2)
            %close()
            
            %disp(size(img_data));
            
            g.putImageQ(img_head, img_data);
            %fprintf('Put on Queue %d, type = %d\n',length(g.Q),g.Q{1}.type);
            
            
            
        end
        
    end
end
