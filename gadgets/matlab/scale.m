classdef scale

    properties

        factor;

    end

    methods

        function obj = config(obj, xmlhdr)
            fprintf(xmlhdr);
        end

        function [reshdr, resdata] = process(obj, head, data)

            reshdr = head;
            reshdr.version = 99;
            resdata = 2*data;

        end

    end
end
