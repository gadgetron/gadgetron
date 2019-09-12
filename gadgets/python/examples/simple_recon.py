
import time
import logging
import ismrmrd

import numpy as np
import numpy.fft as fft

import gadgetron.external

# When running an external Python module, functions decorated with the 'gadgetron.external.module'
# decorator (see the simple_recon function below) will be invoked. This decorated function serves
# as the entry point for external Python modules.


def cfftn(data, axes):
    # Centered fast fourier transform, n-dimensional
    return fft.ifftshift(fft.fftn(fft.fftshift(data, axes=axes), axes=axes, norm='ortho'), axes=axes)


def cifftn(data, axes):
    # Centered inverse fast fourier transform, n-dimensional
    return fft.fftshift(fft.ifftn(fft.ifftshift(data, axes=axes), axes=axes, norm='ortho'), axes=axes)


def noise_adjustment(connection, header):
    # The dataset might include noise measurements (mine does). We'll consume noise measurements, use them
    # to prepare a noise adjustment matrix, and never pass them down the chain. They contain no image data,
    # and will not be missed. We'll also perform noise adjustment on following acquisitions, when we have a
    # noise matrix available.

    noise_matrix = None
    noise_dwell_time = 1.0

    try:
        noise_bandwidth = header.encoding[0].acquisitionSystemInformation.relativeNoiseBandwidth
    except:
        noise_bandwidth = 0.793

    def scaling_factor(acq):
        return np.sqrt(2 * acq.sample_time_us * noise_bandwidth / noise_dwell_time)

    def calculate_whitening_transformation(noise):
        # The details of the whitening transformation is beyond the scope of these comments.
        noise = np.asmatrix(noise)
        covariance = (1.0 / (noise.shape[1] - 1)) * (noise * noise.H)
        return np.linalg.inv(np.linalg.cholesky(covariance))

    def apply_whitening_transformation(acq):
        return np.asarray(scaling_factor(acq) * noise_matrix * np.asmatrix(acq.data))

    def noise_adjust(acq):
        if noise_matrix is not None:
            acq.data[:] = apply_whitening_transformation(acq)
        return acq

    for _, acquisition in connection:
        if acquisition.is_flag_set(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
            noise_matrix = calculate_whitening_transformation(acquisition.data)
            noise_dwell_time = acquisition.sample_time_us
        else:
            yield noise_adjust(acquisition)


def remove_oversampling(iterable, header):
    # The dataset I'm working with was originally taken on a Siemens scanner, and features 2x oversampling
    # along the first image dimension. We're going to have a look at the header. If our encoded space
    # doesn't match the recon space, we're going to crop the acquisition.

    encoding_space = header.encoding[0].encodedSpace.matrixSize
    recon_space = header.encoding[0].reconSpace.matrixSize

    if encoding_space.x == recon_space.x:
        yield from iterable

    x0 = (encoding_space.x - recon_space.x) // 2
    x1 = (encoding_space.x - recon_space.x) // 2 + recon_space.x

    def crop_acquisition(acquisition):
        xspace = cifftn(acquisition.data, axes=[1])
        xspace = xspace[:, x0:x1]
        acquisition.resize(number_of_samples=xspace.shape[1], active_channels=xspace.shape[0])
        acquisition.center_sample = recon_space.x // 2
        acquisition.data[:] = cfftn(xspace, axes=[1])

        return acquisition

    for acq in iterable:
        yield crop_acquisition(acq)


def accumulate_acquisitions(iterable, header):
    # To form images, we need a full slice of data. We accumulate acquisitions until we reach the
    # end of a slice. The acquisitions are then combined in a single buffer, which is then passed
    # on. We also pass on a single acquisition. We don't need the data, but the header is used to
    # initialize some of the image metadata (this is optional, but there's no reason to be cavalier
    # about metadata).

    acquisitions = []
    matrix_size = header.encoding[0].encodedSpace.matrixSize

    def buffer_from_acquisitions(acqs):
        logging.debug(f"Assembling buffer from {len(acqs)} acquisitions.")

        buffer = np.zeros(
            (acqs[0].data.shape[0],
             acqs[0].data.shape[1],
             matrix_size.y,
             matrix_size.z),
            dtype=np.complex64
        )

        for acq in acqs:
            buffer[:, :, acq.idx.kspace_encode_step_1, acq.idx.kspace_encode_step_2] = acq.data

        return buffer

    for acquisition in iterable:
        acquisitions.append(acquisition)
        if acquisition.is_flag_set(ismrmrd.ACQ_LAST_IN_SLICE):
            yield acquisition, buffer_from_acquisitions(acquisitions)
            acquisitions = []


def reconstruct_images(iterable):
    # Reconstruction is an inverse fft. We pass on the reference acquisition unchanged; we'll need
    # it to initialize the image header later.

    for acq, buffer in iterable:
        yield acq, cifftn(buffer, axes=[1, 2, 3])


def combine_channels(iterable):
    # Buffer contains complex images, one for each channel. We combine these into a single image
    # through a sum of squares along the channels (axis 0).

    for acq, buffer in iterable:
        yield acq, np.sqrt(np.sum(np.square(np.abs(buffer)), axis=0))


def create_ismrmrd_images(iterable):
    # Buffer contains the finished image. We wrap it in an ISMRMRD data structure, so the connection
    # can send it back to Gadgetron. We provide a reference acquisition to properly initialize the
    # image header with delicious metadata; feel the good karma.

    for acquisition, buffer in iterable:
        yield ismrmrd.image.Image.from_array(buffer, acquisition=acquisition)


@ gadgetron.external.module
def simple_recon(connection):
    logging.info("Python reconstruction running.")

    start = time.time()

    # Connections are iterable - iterating them can be done once, and will yield the data sent from Gadgetron.
    # In this example, we use a nested sequence of generators (https://wiki.python.org/moin/Generators) each
    # responsible for part of the reconstruction. In this manner, we construct a succession of generators, each
    # one step closer to the final product. Iterating the final iterator thus produces output-ready images.

    iterable = noise_adjustment(connection, connection.header)
    iterable = remove_oversampling(iterable, connection.header)
    iterable = accumulate_acquisitions(iterable, connection.header)
    iterable = reconstruct_images(iterable)
    iterable = combine_channels(iterable)
    iterable = create_ismrmrd_images(iterable)

    # We're only interested in acquisitions in this example, so we filter the connection. Anything filtered out in
    # this way will pass back to Gadgetron unchanged.
    connection.filter(ismrmrd.Acquisition)

    for image in iterable:
        logging.debug("Sending image back to client.")
        connection.send(image)

    logging.info(f"Python reconstruction done. Duration: {(time.time() - start):.2f} s")
