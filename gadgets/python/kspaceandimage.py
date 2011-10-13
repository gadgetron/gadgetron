import numpy as np
import numpy.fft as ft
import numpy.fft.helper as fth

def ktoi(data,axis=-1):
    if (axis == -1):
        ax = fth.arange(0,data.ndim)
    else:
        ax = axis

    return fth.fftshift(ft.ifftn(fth.ifftshift(data,axes=ax),axes=ax),axes=ax)

def itok(data,axis=-1):
    if (axis == -1):
        ax = fth.arange(0,data.ndim)
    else:
        ax = axis


    return fth.fftshift(ft.fftn(fth.ifftshift(data,axes=ax),axes=ax),axes=ax)
