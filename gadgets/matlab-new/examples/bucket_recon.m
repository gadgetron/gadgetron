
function bucket_recon(connection)
    disp("Matlab bucket reconstruction running.") 

    next = slice_from_bucket(@connection.next, connection.header);
    next = reconstruct_slice(next);
    next = combine_channels(next);
    next = create_ismrmrd_image(next);
    next = send_image_to_client(next, connection);
   
    tic, gadgetron.consume(next); toc
end

function next = slice_from_bucket(input, header)

    matrix_size = header.encoding.encodedSpace.matrixSize;

    function image = slice_from_acquisitions(bucket)
        disp("Assembling buffer from bucket containing " + num2str(bucket.data.count) + " acquisitions");
        
        acquisition = ismrmrd.Acquisition(              ...
            single_header(bucket.data.header, 1),       ...
            squeeze(bucket.data.data(1, :, :)),         ...
            squeeze(bucket.data.trajectory(1, :, :))    ...
        );
    
        slice = complex(zeros( ...
            size(acquisition.data, 1), ...
            size(acquisition.data, 2), ...
            matrix_size.y,             ...
            matrix_size.z              ...
        ));
    
        for i = 1:bucket.data.count
            encode_step_1 = bucket.data.header.idx.kspace_encode_step_1(i);
            encode_step_2 = bucket.data.header.idx.kspace_encode_step_2(i);
            slice(:, :, encode_step_1 + 1, encode_step_2 + 1) = squeeze(bucket.data.data(i, :, :));            
        end
        
        disp("Slice sum:")
        s = sum(slice, 'all')
        
        image.data = slice;
        image.acquisition = acquisition;
    end

    next = @() slice_from_acquisitions(input());
end

function header = single_header(header, index)
    take = @(arr) arr(:, index);
    header = structfun(take, header, 'UniformOutput', false);
    header.idx = structfun(take, header.idx, 'UniformOutput', false);
end

% The rest of this file is IDENTICAL to corrosponding functions used in
% simple_recon. They are included here to have self-contained examples, 
% but you should much prefer to have functions like these on the MATLAB
% path and reuse them individually. 

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
