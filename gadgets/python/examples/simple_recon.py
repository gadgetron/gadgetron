
import logging
import ismrmrd

import numpy as np

from multimethod import multimethod, overload

import gadgetron.external

import pyfftw

pyfftw.interfaces.cache.enable()
pyfftw.config.PLANNER_EFFORT = 'FFTW_ESTIMATE'

import pyfftw.interfaces.scipy_fftpack as fftpack


def cfftn(data,axes):
    '''Centered fft'''
    return fftpack.ifftshift(fftpack.fftn(fftpack.fftshift(data,axes=axes)))


def cifftn(data,axes):
    '''Centered ifft'''
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(data,axes),axes=axes))


def calculate_prewhitening(noise, scale_factor=1.0):
    '''Calculates the noise prewhitening matrix
    :param noise: Input noise data (array or matrix), ``[coil, nsamples]``
    :scale_factor: Applied on the noise covariance matrix. Used to
                   adjust for effective noise bandwith and difference in
                   sampling rate between noise calibration and actual measurement:
                   scale_factor = (T_acq_dwell/T_noise_dwell)*NoiseReceiverBandwidthRatio
    :returns w: Prewhitening matrix, ``[coil, coil]``, w*data is prewhitened
    '''

    noise_int = noise.reshape((noise.shape[0], noise.size//noise.shape[0]))
    M = float(noise_int.shape[1])
    dmtx = (1/(M-1))*np.asmatrix(noise_int)*np.asmatrix(noise_int).H
    dmtx = np.linalg.inv(np.linalg.cholesky(dmtx))
    dmtx = dmtx*np.sqrt(2)*np.sqrt(scale_factor)

    return dmtx


def apply_prewhitening(data,dmtx):
    '''Apply the noise prewhitening matrix
    :param noise: Input noise data (array or matrix), ``[coil, ...]``
    :param dmtx: Input noise prewhitening matrix
    :returns w_data: Prewhitened data, ``[coil, ...]``,
    '''

    s = data.shape
    return np.asarray(np.asmatrix(dmtx)*np.asmatrix(data.reshape(data.shape[0],data.size//data.shape[0]))).reshape(s)


def is_noise_adjustment(acq):
    return acq.is_flag_set(ismrmrd.ACQ_IS_NOISE_MEASUREMENT) > 0


def prepare_noise_adjustment(connection):

    noise_matrix = None

    def noise_adjust(acquisition : lambda x: not is_noise_adjustment(x) ):
        if (noise_matrix is not None):
            acquisition.data[:] = apply_prewhitening(acquisition.data,noise_matrix)[:]
        return acquisition

    for mid,acquisition in connection:
        if is_noise_adjustment(acquisition):
            noise_matrix = calculate_prewhitening(acquisition.data)
        else:
            yield noise_adjust(acquisition)


def remove_oversampling(iterable, header: ismrmrd.xsd.ismrmrdHeader):
    # The dataset I'm working with was originally taken on a Siemens scanner, and features 2x oversampling
    # along the first dimension. We're going to have a look at the header. If our encoded space doesn't
    # match the recon space, we're going to crop the acquisitions.

    encoding_space = header.encoding[0].encodedSpace.matrixSize
    recon_space= header.encoding[0].reconSpace.matrixSize

    if encoding_space.x == recon_space.x:
        return iterable

    x0 = (encoding_space.x - recon_space.x) // 2
    x1 = (encoding_space.x - recon_space.x) // 2 + recon_space.x
    for acq in iterable:
        xspace = cifftn(acq.data, axes=[1])
        xspace = xspace[:, x0:x1]
        acq.resize(number_of_samples=xspace.shape[1], active_channels=xspace.shape[0])
        acq.center_sample = recon_space.x // 2
        acq.data[:] = cfftn(xspace, axes=[1])
        yield acq


def accumulate(iterable, header : ismrmrd.xsd.ismrmrdHeader):
    print(header.acquisitionSystemInformation.coilLabel)
    matrixsize = header.encoding[0].encodedSpace.matrixSize
    buffer = None

    is_first = lambda acq : acq.is_flag_set(ismrmrd.ACQ_FIRST_IN_SLICE) > 0
    is_last = lambda acq : acq.is_flag_set(ismrmrd.ACQ_LAST_IN_SLICE) > 0

    def put_acq(acq):
        buffer[:, :, acq.idx.kspace_encode_step_1, acq.idx.kspace_encode_step_2] = acq.data

    @overload
    def handle_acquisition(acq ):
        put_acq(acq)
        return
        yield

    @overload
    def handle_acquisition(acq : is_first):
        nonlocal buffer
        buffer = np.zeros((acq.active_channels,acq.data.shape[1],matrixsize.y,matrixsize.z),dtype=np.complex64)
        put_acq(acq)
        return
        yield

    @overload
    def handle_acquisition(acq : is_last ):
        nonlocal buffer
        put_acq(acq)
        yield buffer

    for acq in iterable:
        yield from handle_acquisition(acq)


def recon(iterable):
    for buffer in iterable:
        yield cifftn(buffer, axes=[1, 2, 3])


def combine_coils(iterable):
    for buffer in iterable:
        rms = np.sqrt(np.sum(np.abs(buffer)**2, axis=0))
        yield rms


def make_image(iterable):
    for buffer in iterable:
        yield ismrmrd.image.Image.from_array(buffer)


@ gadgetron.external.module
def simple_recon(connection):
    logging.info("Python reconstruction running.")

    iterable = prepare_noise_adjustment(connection)
    iterable = remove_oversampling(iterable, connection.header)
    iterable = accumulate(iterable, connection.header)
    iterable = recon(iterable)
    iterable = combine_coils(iterable)
    iterable = make_image(iterable)

    for image in iterable:
        connection.send(image)
