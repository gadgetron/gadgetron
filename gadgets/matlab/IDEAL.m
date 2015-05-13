classdef IDEAL < handle & BaseBufferGadget
    
    properties
        
        image_num;
        series_num;
        center_line;
        img_size;
        frequencies;
        bw;
        dte;
        
    end
    
    methods
        
        function g = config(g)
            
            fprintf('The resonance frequency is %d\n', g.xml.experimentalConditions.H1resonanceFrequency_Hz);
            nx = g.xml.encoding.encodedSpace.matrixSize.x;
            ny = g.xml.encoding.encodedSpace.matrixSize.y;
            nz = g.xml.encoding.encodedSpace.matrixSize.z;
            % for 2D sequences the number of getZ breaks
            %try
            %catch
            
            % nz =1
            %end
            % the number of receiver channels is optional
            try
                % this is the only cast from java.lang.Integer that works in Matlab
                nc = g.xml.acquisitionSystemInformation.receiverChannels;
            catch
                nc = 1;
            end
            % the number of slices is optional
            try
                ns = g.xml.encoding.encodingLimits.slice.maximum + 1;
            catch
                ns = 1;
            end
            
            g.center_line = g.xml.encoding.encodingLimits.kspace_encoding_step_1.center;
            g.img_size = [nx ny nz];
            g.image_num = 0;   % todo this needs to be static or global...
            g.series_num = 0;  % todo this needs to be static or global...
            g.frequencies = [ 392.0   267.0   177.0     0.0  -324.0];
            g.bw = g.xml.userParameters.userParameterDouble(1).value;
            g.dte = g.xml.userParameters.userParameterDouble(2).value;
            %g.frequencies =  [-575.1223 -450.1223 -360.1223 -183.1223  140.8777];
        end
        
        function g = process(g, recon_data)
            disp('Processing')
            % stuff the line
            
            buffer =struct(recon_data(1).data); %We only have the first, and we don't care about reference
            
            [nzf,nte,nslices,ncoils,ntimes]  = size(buffer.data)
            nte = nte-1;
            %% Find frequencies
            dspec =  buffer.data(:,1,:,:,:,:);
            nzf2 = max(nzf,8192);
            spec = fftshift(fft(dspec,nzf2,1),1);
            fax = -g.bw*(-nzf2/2:nzf2/2-1)/nzf2;
            
            
            absspec = sum(sum(sum(sum(abs(spec),2),3),4),5);
            [tmp,imax] = max(absspec);
            %figure;plot(fax,absspec);set(gca,'XDir','reverse');
            %pause(2)
            %pause(2)
            iwidth = ceil(max([2 , 20/(g.bw/nzf2) ])); % width at least +-20Hz or +-2pts
            cog_ind = ((imax-iwidth):(imax+iwidth));
            cog_freq = sum(fax(cog_ind)*absspec(cog_ind,1))/sum(absspec(cog_ind,1));
            fprintf('Pyruvate frequency = %g [Hz]\n',cog_freq);
            freq = g.frequencies+cog_freq
            %freq = -freq
            g.dte
            g.bw
            
            
            t = (0:nte-1).'*g.dte;
            A = exp(-1i*2*pi*t*freq);          % CS kernel
            pinv_A = pinv(A);
            tt = (0:nzf-1)/g.bw;
            
            nfreqs = length(freq);
            %ntimes=1
            outdata = zeros(nzf,nfreqs,nslices,ncoils,ntimes);
            
            
            
            
            for lt=1:ntimes,
                for ls=1:nslices,
                    for lc=1:ncoils
                        
                        
                        outdata(:,:,ls,lc,lt) = reshape(((pinv_A*squeeze(buffer.data(:,2:end,ls,lc,lt)).').* ...
                            exp(1i*2*pi*freq'*tt)).',[nzf,nfreqs]);
                        %                         for lk=1:nzf
                        %                             %buffer.data(lk,2:end,ls,lc,lt)
                        %                             %datatmp(lt,lc,ls,2:end,lk)
                        %
                        %                             outdata(lk,:,ls,lc,lt) = (pinv_A*(buffer.data(lk,2:end,ls,lc,lt))').*exp(-1i*2*pi*tt(lk)*freq');
                        %
                        %
                        %                             %outdata(lk,:,ls,lc,lt) = (pinv_A*reshape(datatmp(lt,lc,ls,2:end,lk),nte,1)).*exp(-1i*2*pi*tt(lk)*freq');
                        %                         end
                    end
                end
            end
            
%             for lf = 1:nfreqs
%                 tmp = outdata(:,lf,:,:,:);
%                 outdata(:,lf,:,:,:) = tmp/norm(tmp(:));
%             end
            
            
            
            %% Create headers
            
            nheaders = nslices*ntimes;
            headers = ismrmrd.AcquisitionHeader(nheaders);
            
            old_headers = buffer.headers;
            headers.version(:) = old_headers.version(1);
            headers.number_of_samples(:) = old_headers.number_of_samples(1);
            headers.center_sample(:) = old_headers.center_sample(1);
            headers.active_channels(:) = ncoils;
            
            headers.read_dir  = repmat(single([1 0 0]'),[1 nheaders]);
            headers.phase_dir = repmat(single([0 1 0]'),[1 nheaders]);
            headers.slice_dir = repmat(single([0 0 1]'),[1 nheaders]);
            
            for rep = 1:ntimes
                acqno = rep;
                
                headers.scan_counter(acqno) = 0;
                headers.idx.kspace_encode_step_1(acqno) = acqno-1;
                headers.idx.repetition(acqno) = 0;
                
                
                % Set the flags
                headers.flagClearAll(acqno);
                if rep == 1
                    
                    headers.flagSet(headers.FLAGS.ACQ_FIRST_IN_ENCODE_STEP1, acqno);
                    headers.flagSet(headers.FLAGS.ACQ_FIRST_IN_SLICE, acqno);
                    headers.flagSet(headers.FLAGS.ACQ_FIRST_IN_REPETITION, acqno);
                end
                if rep==ntimes
                    headers.flagSet(headers.FLAGS.ACQ_LAST_IN_ENCODE_STEP1, acqno);
                    headers.flagSet(headers.FLAGS.ACQ_LAST_IN_SLICE, acqno);
                    headers.flagSet(headers.FLAGS.ACQ_LAST_IN_REPETITION, acqno);
                end
                headers.trajectory_dimensions(acqno) = old_headers.trajectory_dimensions(acqno);
                
                
            end
            
            buffer.headers = headers;
            
            trajsize = size(buffer.trajectory)
            buffer.trajectory = buffer.trajectory(:,:,1,:,:); %Only include the first spirals
            
            
            
            
            %% Create reference data to use for coil sensitivity maps
            reference = struct(buffer); %Hope Matlab does a deep copy
            reference.data = reference.data(:,4,:,:,:); %Only get Pyrovate, which is number 4
            reference.trajectory = reference.trajectory(:,:,1,:,:); %Just get a single spiral
            reference.headers = headers.select(1:ntimes); %Include first header only. Hope it works
            %             reference = recon_data(1).data;
            %             reference.data = reference.data(:,2:end,:,:,:);
            %
            %             reference.trajectory = reference.trajectory(:,:,1:nte);
            % %
            % %             size(buffer.trajectory)
            % %             size(buffer.data)
            %             %% Put on que
            
            for f = 1:nfreqs
                buffer =struct(buffer);
                buffer.data=outdata(:,f,:,:,:);
                g.putBufferQ(buffer,reference);
            end
            %             g.putBufferQ(buffer);
            
            %fprintf('Put on Queue %d, type = %d\n',length(g.Q),g.Q{1}.type);
            
            
            
        end
        
    end
end
