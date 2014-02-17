classdef BaseGadget < handle

    properties

        Q = [];
        xml = [];

    end

    methods

        % Constructor
        function g = BaseGadget()
        end

        % Init function
        function init(g, xmlstr)
            % Convert the xml config string to an IsmrmrdHeader object
            g.xml = org.ismrm.ismrmrd.XMLString.StringToIsmrmrdHeader(xmlstr);
            g.emptyQ();
        end

        % Process function
        function [Q] = run_process(g, htype, hdr_bytes, data)
            if (htype == 1)
                head = ismrmrd.AcquisitionHeader(hdr_bytes);
            elseif (htype == 2)
                head = ismrmrd.ImageHeader(hdr_bytes);
            else
                error('Uknown header type.');
            end
            g.process(head, data);
            Q = g.Q;
        end

        % Config function
        function config(g)
            fprintf('%s\n',char(org.ismrm.ismrmrd.xmlhdr.XMLString.IsmrmrdHeaderToString(g.xml)));
        end
        
        % Process function
        function process(g, head, data)
            g.putQ(head,data);
        end

        % Q related functions
        function emptyQ(g)
           g.Q = [];
        end

        function putQ(g, head, data)
            % find the end of the queue
	        idx = length(g.Q) + 1;
            % put the type of the header and the bytes for the header on the queue
            if isa(head, 'ismrmrd.AcquisitionHeader')
                g.Q(idx).type = int32(1);
                g.Q(idx).bytes = head.toBytes();
            elseif isa(head, 'ismrmrd.ImageHeader')
                g.Q(idx).type = int32(2);
                g.Q(idx).bytes = head.toBytes();
            else
                % TODO: do we throw an error here?
                g.Q(idx).type = int32(0);
            end
            % put the data on the queue
            % make sure the data is single precision
            g.Q(idx).data = single(data);
        end

    end
end
