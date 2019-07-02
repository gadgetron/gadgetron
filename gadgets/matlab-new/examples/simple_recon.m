
function simple_recon(connection)
    disp("Matlab reconstruction running.") 

    tic
    
    next = @connection.next; has_next = @connection.has_next;
    
    [next, has_next] = noise_adjust(next, has_next);
    [next, has_next] = remove_oversampling(next, has_next, connection.header);
    [next, has_next] = accumulate_slice(next, has_next, connection.header);
    [next, has_next] = reconstruct_slice(next, has_next);
    [next, has_next] = combine_channels(next, has_next);
    [next, has_next] = create_ismrmrd_image(next, has_next);
    
    connection.filter('ismrmrd.Acquisition')
    
    while has_next(), connection.send(next()); end
    
    toc
end


function [next, has_next] = noise_adjust(input, has_next)

    noise_matrix = [];

    function transformation = calculate_whitening_transformation(data)
        noise = reshape(data, size(data, 1), []);
        M = single(size(noise, 1));
        transformation = (1/(M-1)) * (noise * noise');
        transformation = inv(chol(transformation));
    end

    function data = apply_whitening_transformation(data)
        if isempty(noise_matrix), return; end
        shape = size(data);
        data = reshape(noise_matrix * reshape(data, shape(1), []), shape);
    end

    function out = handle_noise(acquisition)
        if acquisition.is_flag_set(ismrmrd.Flags.ACQ_IS_NOISE_MEASUREMENT)
            noise_matrix = calculate_whitening_transformation(acquisition.data);
            out = handle_noise(input()); % Recursive call; we output the next item.
        else
            acquisition.data = apply_whitening_transformation(acquisition.data);
            out = acquisition;
        end
    end

    next = @() handle_noise(input());
end

function [next, has_next] = remove_oversampling(input, has_next, header)

    encoding_space = header.encoding.encodedSpace.matrixSize;
    recon_space = header.encoding.reconSpace.matrixSize;

    if encoding_space.x == recon_space.x, next = input; return, end

    x0 = (encoding_space.x - recon_space.x) / 2;
    x1 = (encoding_space.x - recon_space.x) / 2 + recon_space.x;
    along_x_dimension = 2;
    
    function acquisition = remove_oversampling(acquisition)
        xspace = cifft(acquisition.data, along_x_dimension); 
        xspace = xspace(:, x0:x1);
        acquisition.header.number_of_samples = recon_space.x;
        acquisition.header.center_sample = recon_space.x / 2;
        acquisition.data = cfft(xspace, along_x_dimension);
    end

    next = @() remove_oversampling(input());
end

function [next, has_next] = accumulate_slice(input, has_next, header)

    matrix_size = header.encoding.encodedSpace.matrixSize;

    function [slice, acquisition] = slice_from_acquisitions(acquisitions)
        disp("Assembling buffer from " + num2str(length(acquisitions)) + " acquisitions");

        acquisition = head(acquisitions);
       
        slice = complex(zeros( ...
            size(acquisition.data, 1), ...
            size(acquisition.data, 2), ...
            matrix_size.y,             ...
            matrix_size.z              ...
        ));
    
        for acq = acquisitions
            slice(:, :, acq.header.idx.kspace_encode_step_1 + 1, acq.header.idx.kspace_encode_step_2 + 1) = acq.data;
        end
    end

    function slice = accumulate()
        
        acquisitions = gadgetron.util.List.empty;
        
        while has_next()
            acquisition = input();
            acquisitions = cons(acquisitions, acquisition);
            if acquisition.is_flag_set(ismrmrd.Flags.ACQ_LAST_IN_SLICE), break; end
        end
        
        slice = slice_from_acquisitions(acquisitions);
    end

    next = @accumulate;
end

function [next, has_next] = reconstruct_slice(input, has_next)
    
    function [image, acquisition] = reconstruct(slice, acquisition)
        image = cifftn(slice);
    end

    next = @() reconstruct(input());
end

function [next, has_next] = combine_channels(input, has_next)

    function x = square(x), x = x .^ 2; end
    function [image, acquisition] = combine_channels(image, acquisition)
        image = sqrt(sum(square(abs(image)), 1));
    end

    next = @() combine_channels(input());
end

function [next, has_next] = create_ismrmrd_image(input, has_next)
    next = @() input();
end

function data = cfft(data, dim)
    data = ifftshift(fft(fftshift(data, dim), dim), dim);
end

function data = cifft(data, dim)
    data = fftshift(ifft(ifftshift(data, dim), dim), dim);
end

function data = cifftn(data)
    data = fftshift(ifftn(ifftshift(data)));
end
