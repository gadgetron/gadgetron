
function simple_recon(connection)
    disp("Matlab reconstruction running.") 
    
    next = noise_adjust(@connection.next);
    next = remove_oversampling(next, connection.header);
    next = accumulate_slice(next, connection.header);
    next = reconstruct_slice(next);
    next = combine_channels(next);
    next = create_ismrmrd_image(next);
    next = send_image_to_client(next, connection);
    
    connection.filter('ismrmrd.Acquisition')

    tic, gadgetron.consume(next); toc
end


function next = noise_adjust(input)

    noise_matrix = [];

    function transformation = calculate_whitening_transformation(data)
        noise = reshape(data, size(data, 1), []);
        M = size(noise, 2);
        transformation = (1/(M-1)) .* (noise * noise');
        transformation = sqrt(2) .* inv(chol(transformation, 'lower'));
    end

    function data = apply_whitening_transformation(data)
        if isempty(noise_matrix), return; end
        shape = size(data);
        data = reshape(noise_matrix * reshape(data, shape(1), []), shape);
    end

    function acquisition = handle_noise(acquisition)       
        if acquisition.is_flag_set(ismrmrd.Flags.ACQ_IS_NOISE_MEASUREMENT)
            noise_matrix = calculate_whitening_transformation(acquisition.data);
            acquisition = handle_noise(input()); % Recursive call; we output the next item.
        else
            acquisition.data = apply_whitening_transformation(acquisition.data);
        end
    end

    next = @() handle_noise(input());
end

function next = remove_oversampling(input, header)

    encoding_space = header.encoding.encodedSpace.matrixSize;
    recon_space = header.encoding.reconSpace.matrixSize;

    if encoding_space.x == recon_space.x, next = input; return, end

    x0 = (encoding_space.x - recon_space.x) / 2 + 1;
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

function next = accumulate_slice(input, header)

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
    
        for acq = acquisitions.asarray
            slice(:, :, acq.header.idx.kspace_encode_step_1 + 1, acq.header.idx.kspace_encode_step_2 + 1) = acq.data;
        end
    end

    function image = accumulate()
        
        acquisitions = gadgetron.util.List.empty;
        
        while true
            acquisition = input();
            acquisitions = cons(acquisitions, acquisition);
            if acquisition.is_flag_set(ismrmrd.Flags.ACQ_LAST_IN_SLICE), break; end
        end
        
        [image.data, image.acquisition] = slice_from_acquisitions(acquisitions);
    end

    next = @accumulate;
end

function next = reconstruct_slice(input)
    
    function image = reconstruct(image)
        image.data = cifftn(image.data, [2, 3, 4]);
    end

    next = @() reconstruct(input());
end

function next = combine_channels(input)

    function x = square(x), x = x .^ 2; end
    function image = combine_channels(image)
        image.data = sqrt(sum(square(abs(image.data)), 1));
    end

    next = @() combine_channels(input());
end

function next = create_ismrmrd_image(input)

    function image = create_image(image)
        image = ismrmrd.Image.from_data(image.data, image.acquisition);
        image.header.data_type = ismrmrd.Image.DOUBLE;
        image.header.image_type = ismrmrd.Image.MAGNITUDE;
    end

    next = @() create_image(input());
end

function next = send_image_to_client(input, connection)

    function send_image(image)
        disp("Sending image to client.");
        connection.send(image);
    end

    next = @() send_image(input());
end


function data = cfft(data, dim)
    data = ifftshift(fft(fftshift(data, dim), [], dim), dim);
end

function data = cifft(data, dim)
    data = fftshift(ifft(ifftshift(data, dim), [], dim), dim);
end

function data = cifftn(data, dims)
    for dim = dims, data = cifft(data, dim); end
end
