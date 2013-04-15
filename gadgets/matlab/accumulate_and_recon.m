classdef accumulate_and_recon < handle & BaseGadget

    properties

        image_num;
        series_num;
        center_line;
        accumulation;
        
    end

    methods

        function g = config(g)
            fprintf('The resonance frequency is %d\n', g.xml.getExperimentalConditions().getH1ResonanceFrequencyHz());
            nx = g.xml.getEncoding().get(0).getEncodedSpace().getMatrixSize().getX();
            ny = g.xml.getEncoding().get(0).getEncodedSpace().getMatrixSize().getY();
            nz = g.xml.getEncoding().get(0).getEncodedSpace().getMatrixSize().getZ();
            % the number of receiver channels is optional
            try
                % this is the only cast from java.lang.Integer that works in Matlab
                nc = double(g.xml.getAcquisitionSystemInformation().getReceiverChannels());
            catch
	        nc = 1;
            end
            % the number of slices is optional
            try
                ns = g.xml.getEncoding().get(0).getEncodingLimits().getSlice().getMaximum() + 1;
            catch
	        ns = 1;
            end

            g.center_line = g.xml.getEncoding().get(0).getEncodingLimits().getKspaceEncodingStep1().getCenter();
            g.accumulation = zeros(nx, ny, nz, ns, nc);
            g.image_num = 0;  % todo this needs to be static or global...
            g.series_num = 0;  % todo this needs to be static or global...
        end

        function g = process(g, head, data)
            % stuff the line
            line_offset = floor(size(g.accumulation,2)/2) - g.center_line;
            kyind = head.idx.kspace_encode_step_1 + line_offset + 1;
            kzind = head.idx.kspace_encode_step_2 + 1;
            slind = head.idx.slice + 1;
            %fprintf('  offset = %d, center = %d, index = %d\n', line_offset, g.center_line, kyind);

            g.accumulation(:, kyind, kzind, slind, :) = data;

            % At the end of the acquisition, reconstruct the slice
            if (head.isFlagSet('ACQ_LAST_IN_SLICE'))

                img_head = ismrmrd.ImageHeader;
                img_head.channels = head.active_channels;
                img_head.slice = head.idx.slice;
                % set the matrix size
                % note: for 2D size(g.accumulation is only 2D and breaks the setMatrix_size function
                [nx ny nz] = size(g.accumulation);
                img_head.matrix_size = [nx ny nz];
		
                img_head.position = head.position;
                img_head.read_dir = head.read_dir;
                img_head.phase_dir = head.phase_dir;
                img_head.slice_dir = head.slice_dir;
                img_head.patient_table_position = head.patient_table_position;
                img_head.acquisition_time_stamp = head.acquisition_time_stamp;
                img_head.image_index = g.image_num;
                img_head.image_series_index = g.series_num;

		img_data = squeeze(g.accumulation(:,:,:,slind,:));
                img_data = fftshift(ifftn(fftshift(img_data)));
                imagesc(abs(img_data(:,:,1,1))); axis image; axis square;
		pause(2)
                close()

                g.putQ(img_head, img_data);
                %fprintf('Put on Queue %d, type = %d\n',length(g.Q),g.Q{1}.type);

            end

        end

    end
end
