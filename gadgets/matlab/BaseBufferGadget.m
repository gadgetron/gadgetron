classdef BaseBufferGadget < handle
    
  properties

        imageQ = [];
        bufferQ = {};
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
            
            
            for n = 1:numel(recon_data)
                nheaders = numel(recon_data{n}.headers);
                for k = 1:nheaders
                    recon_data{n}.headers{k} = ismrmrd.AcquisitionHeader(recon_data{n}.headers{k});
                end
            end
            
            recon_data{1}.samplingdescription
            g.process(recon_data);
            imageQ = g.imageQ;
            bufferQ = g.bufferQ;
        end
        
          % Q related functions
        function emptyQ(g)
           g.bufferQ = {};
           g.imageQ = [];
        end

        
         function putImageQ(g, header, image)
             disp('Putting Image on Q')
             idx = length(g.imageQ) + 1;
             header.check(); 
             g.imageQ(idx).bytes = header.toBytes();
             g.imageQ(idx).image = single(image);             
         end
        
         function putBufferQ(g,bufferstruct)
             idx = length(g.bufferQ)+1;
             g.bufferQ{idx} = bufferstruct;

             nheaders = length(g.bufferQ{idx}.headers);
             for k = 1:nheaders
                 g.bufferQ{idx}.headers{k} = g.bufferQ{idx}.headers{k}.toBytes();
             end
         end
    end
end