from gadgetron import *
import ismrmrd
import ismrmrd.xsd
import numpy as np
from scipy.fftpack import fft, fftshift
import pylab


class IDEAL(Gadget):
    def __init__(self, next_gadget=None):
        Gadget.__init__(self,next_gadget)
        self.frequencies = []
        self.header = None
        self.bw = None
        self.dte = None


    def process_config(self, conf):
        self.header = ismrmrd.xsd.CreateFromDocument(conf)
        self.enc = self.header.encoding[0]
        print self.header.userParameters
        self.bw = self.header.userParameters.userParameterDouble[0].value_
        self.dte = self.header.userParameters.userParameterDouble[1].value_
        self.frequencies = np.array([float(i) for i in self.get_parameter("frequencies").split(",")])
    def process(self, recon_data):
        buffer = recon_data[0].data
        print np.shape(buffer.data)

        bshape = np.shape(buffer.data)
        bdata = buffer.data
        nzf,nte,nslices,ncoils,ntimes,_,_ = np.shape(bdata) 
        print "input buffer shape",np.shape(bdata)
        nte = nte-1
        dspec = bdata[:,0,...]
        nzf2 = max(nzf,8192)
        spec = fftshift(fft(dspec,n=nzf2,axis=0))
        absspec = np.sum(np.abs(spec),axis=(1,2,3,4))
        imax = np.argmax(absspec)
        iwidth = int(np.floor(max(2,20.0/(self.bw/nzf2))))
        fax = -np.linspace(-0.5,0.5,nzf2)*self.bw
        cog_ind = range(imax-iwidth,imax+iwidth)
        cog_freq = np.sum(fax[cog_ind]*absspec[cog_ind])/np.sum(absspec[cog_ind])
        freqs = self.frequencies+cog_freq

        print "Frequency shift", cog_freq
        t = np.linspace(0,nte*self.dte,nte)
        A = np.matrix(np.exp(-1j*2*np.pi*np.outer(t,freqs)))
        pinv_A = np.linalg.pinv(A)
        tt = np.linspace(0.0,nzf/self.bw,nzf)

        nfreqs = len(freqs)
        outdata = np.zeros((nzf,nfreqs,nslices,ncoils,ntimes),dtype=np.complex64)

        for lt in range(ntimes):
            for ls in range(nslices):
                for lc in range(ncoils):
                    outdata[:,:,ls,lc,lt] = np.reshape(np.array(pinv_A*np.squeeze(bdata[:,1:,ls,lc,lt]).T)*np.exp(1j*2*np.pi*np.outer(freqs,tt)),(nzf,nfreqs))
        nheaders = nslices*ntimes
        baseheader = buffer.headers.flat[0]
        headers = []
        for rep in range(ntimes):
            header = baseheader 
            acqno = rep
            header.scan_counter = 0
            header.idx.kspace_encode_step_1 = acqno
            header.idx.repetition = 0

            if rep == 0:
                header.setFlag(ismrmrd.constants.ACQ_FIRST_IN_ENCODE_STEP1)
                header.setFlag(ismrmrd.constants.ACQ_FIRST_IN_SLICE)
                header.setFlag(ismrmrd.constants.ACQ_FIRST_IN_REPETITION)
            if rep == (ntimes-1):
                header.setFlag(ismrmrd.constants.ACQ_LAST_IN_ENCODE_STEP1)
                header.setFlag(ismrmrd.constants.ACQ_LAST_IN_SLICE)
                header.setFlag(ismrmrd.constants.ACQ_LAST_IN_REPETITION)
            headers.append(header)
        buffer.headers = np.array(headers,dtype=np.dtype(object))


        outdata = outdata[:,:,:,:,:,np.newaxis,np.newaxis]
        ref_data = np.array(outdata[:,3,:,:,:,:,np.newaxis])

        print "Ref data nan",np.sum(ref_data)
        ref_traj = np.array(buffer.trajectory[:,:,0,...])
        ref_headers = buffer.headers.flat[:ntimes]

        buffer.trajectory = buffer.trajectory[:,:,0,0,0,0]
        reference = IsmrmrdDataBuffered(ref_data,ref_headers,sampling=buffer.sampling, trajectory=ref_traj)

        for f in range(nfreqs):
            tmp_buffer = buffer
            tmp_buffer.data = outdata[:,f,np.newaxis,:,:,:,:,:]
            self.put_next([IsmrmrdReconBit(tmp_buffer,reference)])
            #self.put_next([IsmrmrdReconBit(tmp_buffer)])

        return 0





        
            








