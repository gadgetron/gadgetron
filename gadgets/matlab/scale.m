classdef scale < BaseGadget

    properties
        factor;
    end

    methods

        function config(g)
            g.factor = 2;
        end

        function process(g, head, data)
    	    %fprintf('Processing line = %d\n', head.idx.kspace_encode_step_1);
            reshdr = head;
            resdata = g.factor * data;
            g.putQ(reshdr, resdata);
        end

    end
end
