import sys
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
import scipy
from scipy import ndimage

def cpad_2d(data, RO, E1):
    '''
    data: [dRO, dE1, N], padded it round center to [RO, E1, N]
    return value: (s_ro, s_e1), starting point of data in padded array
    '''
    
    try:
        dRO, dE1, N = data.shape
    except:
        data = np.expand_dims(data, axis=2)
        dRO, dE1, N = data.shape

    s_ro = int((RO-dRO)/2)
    s_e1 = int((E1-dE1)/2)
       
    #print(data.shape, RO, E1, s_ro, s_e1)
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
           
    #print(data.shape)
    
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
    """
    Compute landmark from prob
    prob : [RO, E1]
    mode : mean or max
    return pt : [x, y]
    """

    pt = None

    if(binary_mask):
        ind = np.where(prob==thres)
    else:
        if(thres>0 and np.max(prob)<thres):
            return pt
        else:
            mask = adaptive_thresh_cpu(torch.from_numpy(prob), p_thresh=np.max(prob)/2)
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

def get_landmark_from_prob_fast(prob, thres=0.1):
    """
    Compute landmark from prob
    prob : [RO, E1]
    mode : mean or max
    return pt : [x, y]
    """

    #print("get_landmark_from_prob_fast, start ", file=sys.stderr)
    sys.stderr.flush()

    pt = None
    
    if(thres>0 and np.max(prob)<thres):
        return pt
    else:
        v = np.unravel_index(np.argmax(prob), prob.shape)
        pt = np.zeros(2)
        pt[0] = v[1]
        pt[1] = v[0]

        return pt

    #print("get_landmark_from_prob_fast, end ", file=sys.stderr)
    sys.stderr.flush()
    return pt

def adaptive_thresh_cpu(probs, p_thresh=0.5, p_thresh_max=0.988):
    # Try regular adaptive thresholding first
    #p_thresh_max  = 0.988 # <-- Should not be too close to 1 to ensure while loop does not go over.

    p_thresh_incr = 0.01
    #p_thresh = 0.5

    RO = probs.shape[0]
    E1 = probs.shape[1]

    try:
        number_of_blobs = float("inf")
        blobs = np.zeros((RO,E1))
        while number_of_blobs > 1 and p_thresh < p_thresh_max:
            mask = (probs > torch.max(probs) * p_thresh).float()
            blobs, number_of_blobs = ndimage.label(mask)
            p_thresh += p_thresh_incr  # <-- Note this line can lead to float drift.
    
        if(number_of_blobs == 1):
            return mask.numpy()

        if(number_of_blobs == 0):
            mask = np.zeros((RO, E1))
            #print("adaptive_thresh_cpu, did not find any blobs ... ", file=sys.stderr)
            sys.stderr.flush()
            return mask

        ## If we are here then we cannot isolate a singular blob as the LV.
        ## Select the largest blob as the final mask.
        biggest_blob = (0, torch.zeros(RO,E1))
        for i in range(number_of_blobs):
            one_blob = torch.tensor((blobs == i+1).astype(int), dtype=torch.uint8)
            area = torch.sum(one_blob)
            if(area > biggest_blob[0]):
                biggest_blob = (area, one_blob)

        return biggest_blob[1]

    except Exception as e:
        print("Error happened in adaptive_thresh_cpu ...")
        print(e)
        sys.stderr.flush()
        mask = np.zeros((RO,E1))

    return mask
