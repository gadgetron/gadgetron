
function bucket_recon(connection)
    disp("Matlab bucket reconstruction running.") 

    next = create_slice_from_bucket(@connection.next, connection.header);
    next = reconstruct_slice(next);
    next = combine_channels(next);
    next = create_ismrmrd_image(next);
    next = send_image_to_client(next, connection);
   
    tic, gadgetron.consume(next); toc
end

function next = create_slice_from_bucket(input, header)

    matrix_size = header.encoding.encodedSpace.matrixSize;

    function image = slice_from_bucket(bucket)
        disp("Assembling buffer from bucket containing " + num2str(bucket.data.count) + " acquisitions");
        
        image.acquisition = ismrmrd.Acquisition(              ...
            single_header(bucket.data.header, 1),       ...
            squeeze(bucket.data.data(1, :, :)),         ...
            squeeze(bucket.data.trajectory(1, :, :))    ...
        );
    
        image.slice = complex(zeros( ...
            size(image.acquisition.data, 1), ...
            size(image.acquisition.data, 2), ...
            matrix_size.y,             ...
            matrix_size.z              ...
        ));
    
        for i = 1:bucket.data.count
            encode_step_1 = bucket.data.header.idx.kspace_encode_step_1(i);
            encode_step_2 = bucket.data.header.idx.kspace_encode_step_2(i);
            image.slice(:, :, encode_step_1 + 1, encode_step_2 + 1) = ...
                squeeze(bucket.data.data(i, :, :));            
        end
    end

    next = @() slice_from_bucket(input());
end

function header = single_header(header, index)
    take = @(arr) arr(:, index);
    header = structfun(take, header, 'UniformOutput', false);
    header.idx = structfun(take, header.idx, 'UniformOutput', false);
end

% The following functions are IDENTICAL to those found in simple_recon.
% They are repeated here to have self-contained examples, but copy-pasting
% code is of course not adviced. You should keep reusable functions in 
% individual files on your MATLAB path, and reuse them whenever possible.

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

function data = cifft(data, dim)
    data = fftshift(ifft(ifftshift(data, dim), [], dim), dim);
end

function data = cifftn(data, dims)
    for dim = dims, data = cifft(data, dim); end
end
