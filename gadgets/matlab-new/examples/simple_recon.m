
function simple_recon(connection)
    disp("Matlab reconstruction running.") 

    next = @connection.next; has_next = @connection.has_next;
    
    [next, has_next] = noise_adjust(next, has_next);
    [next, has_next] = remove_oversampling(next, has_next, connection.header);
    [next, has_next] = accumulate_slice(next, has_next, connection.header);
    [next, has_next] = reconstruct_slice(next, has_next);
    [next, has_next] = combine_channels(next, has_next);
    [next, has_next] = create_ismrmrd_image(next, has_next);
    
    connection.filter('ismrmrd.Acquisition')
    
    while has_next(), connection.send(next()); end
end


function [next, has_next] = noise_adjust(input, has_next)

    noise_matrix = [];

    function transformation = calculate_whitening_transformation(data)
        noise = reshape(data, size(data, 1), []);
        M = single(size(noise, 1));
        transformation = (1/(M-1)) * (noise * noise');
        transformation = inv(chol(transformation));
        
        disp("Transformation : ")
        transformation
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

    function acquisition = remove_oversampling(acquisition)
        
    end

    next = @() remove_oversampling(input());
end

function [next, has_next] = accumulate_slice(input, has_next, header)

    function [slice, acquisition] = slice_from_acquisitions(acquisitions)

        acquisition = acquisitions{end}
        
        slice = zeros( ...
            size(acquisition.data, 1), ...
            size(acquisition.data, 2), ...
            size(acquisition.data, 3), ...
            length(acquisitions) ...
        );
        
        size(slice)
    end

    function slice = accumulate()
        
        acquisitions = {};
        
        while has_next()
            acquisition = input();
            acquisitions{end + 1} = acquisition;
            if acquisition.is_flag_set(ismrmrd.Flags.ACQ_LAST_IN_SLICE), break; end
        end
        
        slice = slice_from_acquisitions(acquisitions);
    end

    next = @accumulate;
end

function [next, has_next] = reconstruct_slice(input, has_next)
    
    function image = reconstruct(slice)
        image = cifftn(slice, [1, 2, 3]);
    end

    next = @() reconstruct(input());
end

function [next, has_next] = combine_channels(input, has_next)
    next = @() input();
end

function [next, has_next] = create_ismrmrd_image(input, has_next)
    next = @() input();
end


function data = cfftn(data, dim)
    % Centered fast fourier transform, n-dimensional.
    data = ifftshift(fftn(fftshift(data, dim)));
end

function data = cifftn(data, dim)
    % Centered inverse fast fourier transform, n-dimensional.
    data = fftshift(ifftn(ifftshift(data, dim), dim));
end
