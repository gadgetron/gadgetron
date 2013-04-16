classdef mask_image < handle & BaseGadget

    properties
    end

    methods

        function g = config(g)
        end

        function g = process(g, head, data)
            % put the original data on the Q
            g.putQ(head, data);

            % modify the series number
            head.image_series_index = head.image_series_index + 1;

            % zero out a corner of the image
            data(1:end/2,1:end/2,:) = 0;
            
            % put the modified header and image on the Q
            g.putQ(head,data);

        end

    end
end
