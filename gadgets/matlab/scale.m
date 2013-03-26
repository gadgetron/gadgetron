classdef Scale

    properties

        factor;

    end

    methods

        function obj = config(obj, xmlhdr)

        end

        function y = process(obj, head, data, traj)

            y = ismrmrd.Acquisition();
            y.head = head;
            y.data = 2*data;
            if (nargin == 3)
                y.traj = traj;
            end

        end

    end
end
