classdef accumulate_and_recon < handle & BaseGadget

    properties

        image_num;
        series_num;
        center_line;
        accumulation;
        
    end

    methods

        function obj = config(obj, xmlhdr)
            fprintf('The resonance frequency is %d\n', xmlhdr.getExperimentalConditions().getH1ResonanceFrequencyHz());
            nx = xmlhdr.getEncoding().get(0).getEncodedSpace().getMatrixSize().getX();
            ny = xmlhdr.getEncoding().get(0).getEncodedSpace().getMatrixSize().getY();
            nz = xmlhdr.getEncoding().get(0).getEncodedSpace().getMatrixSize().getZ();
            % the number of receiver channels is optional
            try
                % this is the only cast from java.lang.Integer that works in Matlab
                nc = double(xmlhdr.getAcquisitionSystemInformation().getReceiverChannels());
            catch
	        nc = 1;
            end
            % the number of slices is optional
            try
                ns = xmlhdr.getEncoding().get(0).getEncodingLimits().getSlice().getMaximum() + 1;
            catch
	        ns = 1;
            end

            obj.center_line = xmlhdr.getEncoding().get(0).getEncodingLimits().getKspaceEncodingStep1().getCenter();
            obj.accumulation = zeros(nx, ny, nz, ns, nc);
            obj.image_num = 0;  % todo this needs to be static or global...
            obj.series_num = 0;  % todo this needs to be static or global...
        end

        function obj = process(obj, head, data)
            import org.ismrm.ismrmrd.*;

 	    %fprintf('Processing line = %d\n', head.getIdx().getKspace_encode_step_1());

            % stuff the line
            line_offset = floor(size(obj.accumulation,2)/2) - obj.center_line;
            kyind = head.getIdx().getKspace_encode_step_1() + line_offset + 1;
            kzind = head.getIdx().getKspace_encode_step_2() + 1;
            slind = head.getIdx().getSlice() + 1;
            %fprintf('  offset = %d, center = %d, index = %d\n', line_offset, obj.center_line, kyind);

            obj.accumulation(:, kyind, kzind, slind, :) = data;

            % At the end of the acquisition, reconstruct the slice
            lastflag = FlagBit(AcquisitionFlags.ACQ_LAST_IN_SLICE.swigValue());

            if (lastflag.isSet(head.getFlags()))

                img_head = ImageHeader();
                img_head.setChannels(head.getActive_channels());
                img_head.setSlice(head.getIdx().getSlice());
                % set the matrix size
                % note: for 2D size(obj.accumulation is only 2D and breaks the setMatrix_size function
                [nx ny nz] = size(obj.accumulation);
                img_head.setMatrix_size([nx ny nz]);
		
                img_head.setPosition(head.getPosition());
                img_head.setRead_dir(head.getRead_dir());
                img_head.setPhase_dir(head.getPhase_dir());
                img_head.setSlice_dir(head.getSlice_dir());
                img_head.setPatient_table_position(head.getPatient_table_position());
                img_head.setAcquisition_time_stamp(head.getAcquisition_time_stamp());
                img_head.setImage_index(obj.image_num);
                img_head.setImage_series_index(obj.series_num);

		img_data = squeeze(obj.accumulation(:,:,:,slind,:));
	        %subplot(1,2,1); imagesc(abs(img_data(:,:,1,1))); axis image; axis square;
                %subplot(1,2,2); imagesc(abs(fftshift(ifft2(fftshift(img_data(:,:,1,1)))))); axis image; axis square;
		%pause(100)

                obj.putQ(img_head, img_data);
                %fprintf('Put on Queue %d, type = %d\n',length(obj.Q),obj.Q{1}.type);

            end

        end

    end
end
