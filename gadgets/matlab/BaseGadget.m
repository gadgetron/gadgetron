classdef BaseGadget

    properties

        Q;
        head;
        data;

    end

    methods

        function obj = BaseGadget()
            obj.Q = {};
        end

        function obj = config(obj, xmlhdr)

        end

        function obj = putIn(obj, type, bytes, data)
            import org.ismrm.ismrmrd.*;
            import org.ismrm.ismrmrd.xmlhdr.*;
            fprintf('hello from putIn');
            if type == 1
                head = AcquisitionHeader();
                ismrmrd.copyJBytesToAcquisitionHeader(bytes, head);
            elseif type == 2
                head = ImageHeader();
                ismrmrd.copyJBytesToImageHeader(bytes, head);
            else
                % TODO: throw error
                return;
            end
            obj.head = head;
            obj.data = data;
         end

        function obj = putQ(obj, head, data)
            import org.ismrm.ismrmrd.*;
            import org.ismrm.ismrmrd.xmlhdr.*;
            if isa(head, 'AcquisitionHeader')
                obj.Q{end+1}.type = int32(1);
                obj.Q{end+1}.bytes = ismrmrd.acquisitionHeaderToJBytes(head);
            elseif isa(head, 'ImageHeader')
                obj.Q{end+1}.type = int32(2);
                obj.Q{end+1}.bytes = ismrmrd.imageHeaderToJBytes(head);
            else
                % TODO: throw error
                return;
            end
            obj.Q{end+1}.data = data;
        end

    end
end
