classdef scale < BaseGadget

    properties

        factor;

    end

    methods

        function obj = scale()
            obj = obj@BaseGadget();
            obj.factor = 2;
        end

        function obj = config(obj, xmlhdr)
            import org.ismrm.ismrmrd.*;
            import org.ismrm.ismrmrd.xmlhdr.*;
            ihdr = XMLString.StringToIsmrmrdHeader(xmlhdr);
            fprintf('The resonance frequency is %d\n', ihdr.getExperimentalConditions().getH1ResonanceFrequencyHz());
        end

        function obj = process(obj)

            import org.ismrm.ismrmrd.*;
            import org.ismrm.ismrmrd.xmlhdr.*;
            fprintf(class(obj.head))
            reshdr = obj.head;
            reshdr.setVersion(99);
            resdata = obj.factor * obj.data;

            obj.putQ(reshdr, resdata);

        end

    end
end
