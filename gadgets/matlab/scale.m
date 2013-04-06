classdef scale < BaseGadget

    properties

        factor;
        
    end

    methods

        function obj = config(obj, xmlhdr)
            fprintf('The resonance frequency is %d\n', xmlhdr.getExperimentalConditions().getH1ResonanceFrequencyHz());
            obj.factor = 2;
        end

        function obj = process(obj, head, data)

	    %fprintf('Processing line = %d\n', head.getIdx().getKspace_encode_step_1());
            reshdr = head;
            reshdr.setVersion(99);
            
            resdata = obj.factor * data;

            obj = obj.putQ(reshdr, resdata);
            %fprintf('Put on Queue %d, type = %d\n',length(obj.Q),obj.Q{1}.type);

        end

    end
end
