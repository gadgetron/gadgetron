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
            % this is the only cast from java.lang.Integer that works in Matlab
            nc = double(xmlhdr.getAcquisitionSystemInformation().getReceiverChannels()); 
            ns = xmlhdr.getEncoding().get(0).getEncodingLimits().getSlice().getMaximum() + 1;
            obj.center_line = xmlhdr.getEncoding().get(0).getEncodingLimits().getKspaceEncodingStep1().getCenter();
            obj.accumulation = zeros(nx, ny, nz, ns, nc);
            obj.image_num = 0;  % todo this needs to be static or global...
            obj.series_num = 0;  % todo this needs to be static or global...
        end

        function obj = process(obj, head, data)
            import org.ismrm.ismrmrd.*;
            % stuff the line
            line_offset = floor(size(obj.accumulation,2)/2) - obj.center_line;
            obj.accumulation(:,                                                 ...
                    head.getIdx().getKspace_encode_step_1() + 1 + line_offset,  ...
                    head.getIdx().getKspace_encode_step_2() + 1,                ...
                    head.getIdx().getSlice()+1,                                 ...
                    :) = data;
	    %fprintf('Processing line = %d\n', head.getIdx().getKspace_encode_step_1());

            % At the end of the acquisition, reconstruct the slice
            lastflag = FlagBit(AcquisitionFlags.ACQ_LAST_IN_SLICE.swigValue());

            if (lastflag.isSet(head.getFlags()))
                img_head = ImageHeader();
                img_head.setChannels(head.getActive_channels());
                img_head.setSlice(head.getIdx().getSlice());
                dims = size(obj.accumulation);
                img_head.setMatrix_size(dims(1:3));
                img_head.setPosition(head.getPosition());
                img_head.setRead_dir(head.getRead_dir());
                img_head.setPhase_dir(head.getPhase_dir());
                img_head.setSlice_dir(head.getSlice_dir());
                img_head.setPatient_table_position(head.getPatient_table_position());
                img_head.setAcquisition_time_stamp(head.getAcquisition_time_stamp());
                img_head.setImage_index(obj.image_num);
                img_head.setImage_series_index(obj.series_num);

                obj.putQ(img_head, obj.accumulation);
                %fprintf('Put on Queue %d, type = %d\n',length(obj.Q),obj.Q{1}.type);
            end

        end

    end
end
