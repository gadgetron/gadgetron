classdef BaseGadget < handle

    properties

        Q = {};

    end

    methods

        function obj = config(obj, xmlhdr)

        end

        function obj = process(obj, head, data)

        end

        function obj = emptyQ(obj)
           obj.Q = {};
        end

        function ret = getQLength(obj)
	  ret = int32(length(obj.Q));
        end

        function ret = getQ(obj,idx)
	  ret = obj.Q{idx};
        end

        function obj = putQ(obj, head, data)
            import org.ismrm.ismrmrd.*;
            import org.ismrm.ismrmrd.xmlhdr.*;
            
            idx = length(obj.Q) + 1;
            if isa(head, 'AcquisitionHeader')
	        obj.Q{idx}.type = int32(1);
                obj.Q{idx}.bytes = ismrmrd.acquisitionHeaderToJBytes(head);
            elseif isa(head, 'ImageHeader')
	        obj.Q{idx}.type = int32(2);
                obj.Q{idx}.bytes = ismrmrd.imageHeaderToJBytes(head);
            else
                % TODO: throw error
		obj.Q{idx}.type = int32(0);
                return;
            end
            obj.Q{idx}.data = data;
            
        end

    end
end
