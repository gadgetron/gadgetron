classdef BaseBufferGadget < handle
    
  properties

        imageQ = [];
        bufferQ = struct([]);
        xml = [];

    end

    methods

        % Constructor
        function g = BaseBufferGadget()
            
        end

        % Init function
        function init(g, xmlstr)
            % Convert the xml config string to an IsmrmrdHeader object
            g.xml = ismrmrd.xml.deserialize(xmlstr);
            g.emptyQ();
        end
        
           % Process function
        function [imageQ,bufferQ] = run_process(g, recon_data)
            %% Convert headers
            for n = 1:numel(recon_data)
                recon_data(n).data.headers = ismrmrd.AcquisitionHeader(recon_data(n).data.headers);
                if isfield(recon_data(n),'reference')
                    if isstruct(recon_data(n).reference)
                        recon_data(n).reference.headers = ismrmrd.AcquisitionHeader(recon_data(n).reference.headers);
                    end
                end
            end
            
            %% Process data
            g.process(recon_data);
            
            imageQ = g.imageQ;
            bufferQ = g.bufferQ;
        end
        
          % Q related functions
        function emptyQ(g)
           g.bufferQ = struct([]);
           g.imageQ = [];
        end

        
         function putImageQ(g, header, image)
             disp('Putting Image on Q')
             idx = length(g.imageQ) + 1;
             header.check(); 
             g.imageQ(idx).bytes = header.toBytes();
             g.imageQ(idx).image = single(image);             
         end
        
         function putBufferQ(g,data,reference)
             idx = length(g.bufferQ)+1;
             g.bufferQ(idx).data = data;
             g.bufferQ(idx).data.headers = g.bufferQ(idx).data.headers.toBytes();
             if (exist('reference','var'))
                 g.bufferQ(idx).reference = reference;
                 g.bufferQ(idx).reference.headers = reference.headers.toBytes();
             end
        end
    end
end