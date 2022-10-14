"""Utility file for cmr landmark detection"""

import os
import sys
import urllib.request
import hashlib
import fcntl
import numpy as np
from scipy import ndimage

class model_file_lock:
    """A python inter-process file lock
    
    Adaptive from https://stackoverflow.com/questions/6931342/system-wide-mutex-in-python-on-linux
    """
    
    def __init__(self, model_dest, model_sha256):
        self.model_dest = model_dest
        self.model_sha256 = model_sha256
    
    def __enter__ (self):
        self.fp = open(os.path.join(self.model_dest, self.model_sha256 + "_ai_model_lockfile.lck"), 'w')
        fcntl.flock(self.fp.fileno(), fcntl.LOCK_EX)

    def __exit__ (self, _type, value, tb):
        fcntl.flock(self.fp.fileno(), fcntl.LOCK_UN)
        self.fp.close()
            
def is_model_valid(model_file, model_sha256):
    """Check whether the sha256 hash is valid"""
    
    if not os.path.isfile(model_file):
        return False

    with open(model_file, "rb") as f:
        file_hash = hashlib.sha256()
        while chunk := f.read(65536):
            file_hash.update(chunk)
        
    return model_sha256 == file_hash.hexdigest()


def check_and_get_model(model_host, model_file, model_dest, model_sha256):
    """Get the model from the storage account.
    
    Input:
        model_host [str] : host address to download the model
        model_file [str] : model file name
        model_dest [str] : destination path to store the model
        model_sha256 [str] : model sha256 hash
        
    If the model already exist, check its integrity; if invalid, download file
    If model cannot be downloaded, an exception is thrown
    """
    try:
        print(f"check_and_get_model,  model_host is {model_host}", file=sys.stderr)
        print(f"check_and_get_model,  model_file is {model_file}", file=sys.stderr)
        print(f"check_and_get_model,  model_dest is {model_dest}", file=sys.stderr)
        print(f"check_and_get_model,  model_sha256 is {model_sha256}", file=sys.stderr)
        
        # assemble the source and destination of model
        model_url = model_host + model_file
        destination = os.path.join(model_dest, model_file)
        
        print(f"check_and_get_model,  model_url is {model_url}", file=sys.stderr)
        print(f"check_and_get_model,  destination is {destination}", file=sys.stderr)
                   
        with model_file_lock(model_dest, model_sha256):
            # download the model if not exist
            if(not os.path.isfile(destination)):            
                print(f"check_and_get_model, start downloading model", file=sys.stderr)
                urllib.request.urlretrieve(model_url, destination)
                print(f"check_and_get_model, finish downloading model", file=sys.stderr)

            if not is_model_valid(destination, model_sha256):
                print(f"Downloaded model file {destination} failed in sha256 validation.", file=sys.stderr)
                urllib.request.urlretrieve(model_url, destination)
                    
                if not is_model_valid(destination, model_sha256):
                    print(f"Newly downloaded model file {destination} failed in sha256 validation.", file=sys.stderr)
                    raise "invalid model"
            else:
                print(f"Downloaded model file {destination} succeeded in sha256 validation.", file=sys.stderr)
                
    except Exception as e:
            print("Error happened in check_and_get_model ... ", file=sys.stderr)
            print(e, file=sys.stderr)
            raise e

def cpad_2d(data, RO, E1):
    """center crop/pad image to required size
    
    Inputs:
        data: numpy array [dRO, dE1, N]
        
    Return:
        data_padded: numpy array [RO, E1, N], padded/cropped around the center
        s_ro, s_e1: starting indexes of original array in the padded array
    """
    
    if(data.ndim==2):
        data = np.expand_dims(data, axis=2)       
    dRO, dE1, N = data.shape

    s_ro = int((RO-dRO)/2)
    s_e1 = int((E1-dE1)/2)
       
    if(s_ro>=0):
        data_padded = np.zeros((RO, dE1, N))
        if(dRO>=RO):
            data_padded = data[s_ro:s_ro+RO,:,:]
        else:
            data_padded[s_ro:s_ro+dRO,:,:] = data
        data = data_padded
    else:
        data_padded = np.zeros((RO, dE1, N))
        if(dRO+s_ro+s_ro>RO):
            data_padded = data[-s_ro:(dRO+s_ro-1),:,:]
        else:
            data_padded = data[-s_ro:(dRO+s_ro),:,:]
        data = data_padded
           
    if(s_e1>=0):
        data_padded = np.zeros((RO, E1, N))
        if(dE1>=E1):
            data_padded = data[:,s_e1:s_e1+E1,:]
        else:
            data_padded[:,s_e1:s_e1+dE1,:] = data
        data = data_padded
    else:
        data_padded = np.zeros((RO, E1, N))
        if(dE1+s_e1+s_e1>E1):
            data_padded = data[:,-s_e1:(dE1+s_e1-1),:]
        else:
            data_padded = data[:,-s_e1:(dE1+s_e1),:]
        data = data_padded
       
    return data_padded, s_ro, s_e1

def get_landmark_from_prob(prob, thres=0.5, mode="mean", binary_mask=False):
    """Compute landmark location from the model probablity maps
    
    Inputs:
        prob : [RO, E1], the model produced probablity map for a landmark
        thres : if np.max(prob)<thres, determine there is no landmark detected
        mode : mean or max, use mean or max probablity to detect landmark
        binary_mask : if true, prob is a binary (0 or 1) map
        
    Outputs:
        pt : [x, y], detected landmark point
    """

    pt = None

    if(binary_mask):
        ind = np.where(prob==thres)
    else:
        if(thres>0 and np.max(prob)<thres):
            return pt
        else:
            adaptive_thres = 0.5
            mask = adaptive_thresh_cpu(prob, p_thresh=adaptive_thres*np.max(prob))
            ind = np.where(mask>0)

    if (np.size(ind[0])==0):
        return pt

    pt = np.zeros(2)
    if(mode == "mean"):
        pt[0] = np.mean(ind[1].astype(np.float32))
        pt[1] = np.mean(ind[0].astype(np.float32))
    else:
        v = np.unravel_index(np.argmax(prob), prob.shape)
        pt[0] = v[1]
        pt[1] = v[0]

    return pt

def adaptive_thresh_cpu(probs, p_thresh=0.5, p_thresh_max=0.988):
    """Adaptive thresholding
    
    Gradually increasing the threshold from p_thresh to p_thresh_max, until
    there is only one blob is found in probs
    
    Inputs:
        probs: [RO, E1], probalibity map computed by model
        p_thresh : staring threshold
        p_thresh_max: max threshold
        
    Outputs:    
        mask: [RO, E1], the binary mask for the detected blob
    """

    p_thresh_incr = 0.01

    RO = probs.shape[0]
    E1 = probs.shape[1]

    try:
        number_of_blobs = float("inf")
        blobs = np.zeros((RO,E1))
        while number_of_blobs > 1 and p_thresh < p_thresh_max:
            mask = (probs > np.max(probs) * p_thresh).astype(np.float32)
            blobs, number_of_blobs = ndimage.label(mask)
            p_thresh += p_thresh_incr
    
        if(number_of_blobs == 1):
            return mask

        if(number_of_blobs == 0):
            mask = np.zeros((RO, E1))
            file=sys.stderr.flush()
            return mask

        # if there are more than one blob, try to isolate the high probablity one
        biggest_blob = (0, np.zeros((RO,E1)))
        for i in range(number_of_blobs):
            one_blob = np.zeros_like(blobs)
            one_blob[blobs == i+1] = 1
            area = np.sum(one_blob)
            if(area > biggest_blob[0]):
                biggest_blob = (area, one_blob)

        return biggest_blob[1]

    except Exception as e:
        print("Error happened in adaptive_thresh_cpu ...")
        print(e)
        file=sys.stderr.flush()
        mask = np.zeros((RO,E1))

    return mask

if __name__=="__main__":
    
    GT_HOME = os.environ['GADGETRON_HOME']
    GT_PYTHON_DIR = os.path.join(GT_HOME, 'share/gadgetron/python')
    print("GT_HOME is", GT_HOME)
    
    # test get_model
    model_host = 'https://gadgetrondata.blob.core.windows.net/cmr-ai-models/'
    model_file = 'CMR_landmark_network_RO_320_E1_320_CH2_CH3_CH4_Myo_Pts_sFOV_LossMultiSoftProb_KLD_Dice_Pytorch_1.8.0a0+37c1f4a_2021-08-08_20210808_085042.onnx'
    model_dest = os.path.join(GT_PYTHON_DIR, 'cmr_lax_landmark_detection')
    model_sha256 = '48efe3e70b1ff083c9dd0066469f62bf495e52857d68893296e7375b69f691e4'
    
    res = check_and_get_model(model_host, model_file, model_dest, model_sha256)
