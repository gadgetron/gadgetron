classdef trajectoryScale < handle & BaseBufferGadget
    
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
            
            
           for n = 1:numel(recon_data) 
               buffer =recon_data{n};
               buffer.trajectory = buffer.trajectory*2;
               fieldnames(buffer)
               g.putBufferQ(buffer) 
           end
           
           
           
            
            %fprintf('Put on Queue %d, type = %d\n',length(g.Q),g.Q{1}.type);
            
            
            
        end
        
    end
end
