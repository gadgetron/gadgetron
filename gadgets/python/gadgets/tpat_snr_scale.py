import ismrmrd
import ismrmrd.xsd
import numpy as np
from ismrmrdtools import transform, coils, grappa
from gadgetron import Gadget
import copy 
import math

class RemOS(Gadget):
    def process_config(self, conf):
        return

    def process(self, acq, data,*args):
        if not acq.isFlagSet(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
            ro_length = acq.number_of_samples
            padded_ro_length = (acq.number_of_samples-acq.center_sample)*2
            if padded_ro_length != ro_length: #partial fourier
                data2 = np.zeros((data.shape[0], padded_ro_length),dtype=np.complex64)
                offset = (padded_ro_length>>1)  - acq.center_sample
                data2[:,0+offset:offset+ro_length] = data
            else:
                data2 = data
    
            data2=transform.transform_kspace_to_image(data2,dim=(1,))
            data2=data2[:,(padded_ro_length>>2):(padded_ro_length>>2)+(padded_ro_length>>1)]
            data2=transform.transform_image_to_kspace(data2,dim=(1,)) * np.sqrt(float(padded_ro_length)/ro_length)
            acq.center_sample = padded_ro_length>>2
            acq.number_of_samples = data2.shape[1]
            self.put_next(acq,data2,*args)
        return 0                                                                                     

class NoiseAdj(Gadget):
    def __init__(self, next_gadget = None):
        Gadget.__init__(self, next_gadget)
        self.noise_data = list()
        self.noise_dmtx = None
    def process(self,acq,data,*args):
        if acq.isFlagSet(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
            self.noise_data.append((acq,data))
        else:
            if len(self.noise_data):
                profiles = len(self.noise_data)
                channels = self.noise_data[0][1].shape[0]
                samples_per_profile = self.noise_data[0][1].shape[1]
                noise = np.zeros((channels,profiles*samples_per_profile),dtype=np.complex64)
                counter = 0
                for p in self.noise_data:
                    noise[:,counter*samples_per_profile:(counter*samples_per_profile+samples_per_profile)] = p[1]
                    counter = counter + 1
                
                scale = (acq.sample_time_us/self.noise_data[0][0].sample_time_us)*0.79
                self.noise_dmtx = coils.calculate_prewhitening(noise,scale_factor=scale)
                
                #Test the noise adjust
                d = self.noise_data[0][1]
                d2 = coils.apply_prewhitening(d, self.noise_dmtx)                
                self.noise_data = list()
            
            if self.noise_dmtx is not None:
                data2 = coils.apply_prewhitening(data, self.noise_dmtx)
            else:
                data2 = data
                
            self.put_next(acq,data2)
        return 0

class PCA(Gadget):
    def __init__(self, next_gadget=None):
        Gadget.__init__(self, next_gadget) 
        self.calib_data = list()
        self.pca_mtx = None
        self.max_calib_profiles = 100
        self.samples_to_use = 16
        self.buffering = True
        
    def process(self,acq,data,*args):    
        if self.buffering:
            self.calib_data.append((acq,data))
            
            if (len(self.calib_data)>=self.max_calib_profiles or acq.isFlagSet(ismrmrd.ACQ_LAST_IN_SLICE)):
                #We are done buffering calculate pca transformation
                if self.samples_to_use < acq.number_of_samples:
                    samp_to_use = self.samples_to_use
                    
                if (len(self.calib_data) < 16):
                    samp_to_use = acq.number_of_samples
                    
                total_samples = samp_to_use*len(self.calib_data)
                channels = data.shape[0]
                
                A = np.zeros((total_samples,channels), dtype=np.complex64)
                counter = 0
                for p in self.calib_data:
                    d = p[1][:, acq.center_sample-(samp_to_use>>1):acq.center_sample+(samp_to_use>>1)]
                    A[counter*samp_to_use:counter*samp_to_use+samp_to_use,:] = np.transpose(d)
                    counter = counter+1
                
                m = np.mean(A,0)
                A_m = A - m.reshape((1,m.shape[0]))
                U, s, V = np.linalg.svd(A_m, full_matrices=False)
                
                self.pca_mtx = V
                
                #Empty calib_data
                for p in self.calib_data:
                    data2 = np.dot(self.pca_mtx,p[1])
                    self.put_next(p[0],data2)
     
                self.buffering = False
                self.calib_data = list()
                return 0
        else:
            if self.pca_mtx is not None:
                data2 = np.dot(self.pca_mtx,data)
                self.put_next(acq,data2,*args)
            else:
                self.put_next(acq,data,*args)
            
        return 0

class CoilReduce(Gadget):
    def __init__(self, next_gadget = None):
        Gadget.__init__(self, next_gadget)
        self.coils_out = 16
        
    def process(self, acq, data, *args):
        if acq.active_channels > self.coils_out:
            data2 = data[0:self.coils_out,:]
            acq.active_channels = self.coils_out
        else:
            data2 = data
            
        self.put_next(acq,data2,*args)
        return 0


class Recon(Gadget):    
    def __init__(self, next_gadget=None):
        Gadget.__init__(self, next_gadget) 
        self.header = None
        self.enc = None
        self.acc_factor = None
        self.buffer = None
        self.samp_mask = None
        self.header_proto = None
        self.calib_buffer = list()
        self.unmix = None
        self.gmap = None
        self.calib_frames = 0
    
    def process_config(self, cfg):
        self.header = ismrmrd.xsd.CreateFromDocument(cfg)
        self.enc = self.header.encoding[0]

        #Parallel imaging factor
        self.acc_factor = self.enc.parallelImaging.accelerationFactor.kspace_encoding_step_1
        
        reps = self.enc.encodingLimits.repetition.maximum+1
        phs = self.enc.encodingLimits.phase.maximum+1
        if reps > phs:
            self.calib_frames = reps
        else:
            self.calib_frames = phs
            
        if self.calib_frames < self.acc_factor:
            self.calib_frames = self.acc_factor
        
        #Frames should be a multiple of the acceleration factor
        self.frames = math.floor(self.calib_frames/self.acc_factor)*self.acc_factor
        
    def process(self, acq, data,*args):

        if self.buffer is None:
            # Matrix size
            eNx = self.enc.encodedSpace.matrixSize.x
            eNy = self.enc.encodedSpace.matrixSize.y
            eNz = self.enc.encodedSpace.matrixSize.z
            rNx = self.enc.reconSpace.matrixSize.x
            rNy = self.enc.reconSpace.matrixSize.y
            rNz = self.enc.reconSpace.matrixSize.z

            # Field of View
            eFOVx = self.enc.encodedSpace.fieldOfView_mm.x
            eFOVy = self.enc.encodedSpace.fieldOfView_mm.y
            eFOVz = self.enc.encodedSpace.fieldOfView_mm.z
            rFOVx = self.enc.reconSpace.fieldOfView_mm.x
            rFOVy = self.enc.reconSpace.fieldOfView_mm.y
            rFOVz = self.enc.reconSpace.fieldOfView_mm.z
        
            channels = acq.active_channels

            if data.shape[1] != rNx:
                raise("Error, Recon gadget expects data to be on correct matrix size in RO direction")
                
            if (rNz != 1):
                rasie("Error Recon Gadget only supports 2D for now")
                
            self.buffer = np.zeros((channels, rNy, rNx),dtype=np.complex64)
            self.samp_mask = np.zeros(self.buffer.shape[1:])
            self.header_proto = ismrmrd.ImageHeader()
            self.header_proto.matrix_size[0] = rNx
            self.header_proto.matrix_size[1] = rNy
            self.header_proto.matrix_size[2] = rNz
            self.header_proto.field_of_view[0] = rFOVx
            self.header_proto.field_of_view[1] = rFOVy
            self.header_proto.field_of_view[0] = rFOVz
        
        #Now put data in buffer
        line_offset = self.buffer.shape[1]/2 - self.enc.encodingLimits.kspace_encoding_step_1.center                                                                                 
        self.buffer[:,acq.idx.kspace_encode_step_1+line_offset,:] = data                                                          
        self.samp_mask[acq.idx.kspace_encode_step_1+line_offset,:] = 1
        
        #If last scan in buffer, do FFT and fill image header
        if acq.isFlagSet(ismrmrd.ACQ_LAST_IN_ENCODE_STEP1) or acq.isFlagSet(ismrmrd.ACQ_LAST_IN_SLICE):
            img_head = copy.deepcopy(self.header_proto)
            img_head.position = acq.position                                                                                                                               
            img_head.read_dir = acq.read_dir                                                                                                                               
            img_head.phase_dir = acq.phase_dir                                                                                                                             
            img_head.slice_dir = acq.slice_dir                                                                                                                             
            img_head.patient_table_position = acq.patient_table_position                                                                                                   
            img_head.acquisition_time_stamp = acq.acquisition_time_stamp                                                                                                   
            img_head.slice = acq.idx.slice
            img_head.channels = 1
            
            scale = self.samp_mask.size/(1.0*np.sum(self.samp_mask[:]));

            #We have not yet calculated unmixing coefficients
            if self.unmix is None:
                self.calib_buffer.append((img_head,self.buffer.copy()))
                self.buffer[:] = 0
                self.samp_mask[:] = 0
                
                if len(self.calib_buffer) >= self.calib_frames:
                    cal_data = np.zeros(self.calib_buffer[0][1].shape, dtype=np.complex64)
                    for c in self.calib_buffer:
                        cal_data = cal_data + c[1]
                        
                    mask = np.squeeze(np.sum(np.abs(cal_data),0))
                    mask = np.ones(mask.shape)*(np.abs(mask)>0.0)
                    target = None #cal_data[0:8,:,:]
                    self.unmix, self.gmap = grappa.calculate_grappa_unmixing(cal_data, self.acc_factor, data_mask=mask, kernel_size=(4,5), target_data=target)

                    for c in self.calib_buffer:
                        recon = transform.transform_kspace_to_image(c[1],dim=(1,2))*np.sqrt(scale)
                        recon = np.squeeze(np.sum(recon * self.unmix,0))
                        self.put_next(c[0], recon,*args)
                        
                return 0
                
            if self.unmix is None:
                raise Exception("We should never reach this point without unmixing coefficients")
                
            recon = transform.transform_kspace_to_image(self.buffer,dim=(1,2))*np.sqrt(scale)
            recon = np.squeeze(np.sum(recon * self.unmix,0))
            self.buffer[:] = 0
            self.samp_mask[:] = 0
            self.put_next(img_head,recon,*args)
        return 0
