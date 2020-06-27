from __future__ import print_function

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torchvision.transforms as T
from torchvision.utils import *
import numpy as np
import scipy
import time
import imp
import os
import sys
import math
import copy
import random
import shutil
import scipy.misc
from scipy.io import loadmat, savemat
from glob import glob
import logging

def complex_2_real(a):
    a_r = np.real(a)
    a_i = np.imag(a)

    D, NA = a_r.shape

    res = np.zeros((D, 2*NA), dtype=np.float32)

    res[:,0:NA] = a_r
    res[:,NA:] = a_i

    return res

# -------------------------------------------------------------------------
def real_2_complex(a):
    D, NA = a.shape

    N = int(NA/2)
    res = a[:,0:N] + 1j * a[:,N:]
    return res

# -------------------------------------------------------------------------

class GrappaAI(nn.Module):
    def __init__(self, Din, Dout, params):
        
        super(GrappaAI, self).__init__()

        self.verbose = params['verbose']
        self.with_bias = False
        self.grappa_weights_r = params['grappa_weights_r']
        self.grappa_weights_i = params['grappa_weights_i']
        
        # for real and imag
        self.input_layer_r = torch.nn.Linear(int(Din/2), int(Dout/2), bias=self.with_bias)
        self.input_layer_i = torch.nn.Linear(int(Din/2), int(Dout/2), bias=self.with_bias)

        if self.verbose:
            print("    GrappaAI : input size (%d), output size (%d), with bias %s" % (Din, Dout, self.with_bias))

        print("---> Set grappa weights")
        self.input_layer_r.weight.data = self.grappa_weights_r
        self.input_layer_i.weight.data = self.grappa_weights_i
                
    def forward(self, x):

        # split x to real and img
        M, N = x.shape

        N_h = int(N/2)

        x_r = x[:, 0:N_h]
        x_i = x[:, N_h:N]

        r_r = self.input_layer_r(x_r)
        i_i = self.input_layer_i(x_i)

        r_i = self.input_layer_r(x_i)
        i_r = self.input_layer_i(x_r)

        res_r = r_r - i_i
        res_i = r_i + i_r

        out = torch.cat([res_r, res_i], dim=1)

        return out

# -------------------------------------------------------------------------

def create_grappa_ai_model(gt_ker):
    '''gt_ker is the grappa calibration kernel of A*ker = B

    Dout is the number of output elements. It is 2*B.shape(1), due to complex numbers
    '''

    model = None

    try:
        print("Grappa ai", file=sys.stderr)
        print("----------------------------------------", file=sys.stderr)

        try:
            kRO, kNE1, srcCHA, dstCHA, oE1 = gt_ker.shape
            gt_ker_2D = np.reshape(gt_ker, (kRO*kNE1*srcCHA, dstCHA*oE1), order='F')
        except:
            gt_ker_2D = gt_ker

        gt_ker_2D_r = np.real(gt_ker_2D)
        gt_ker_2D_i = np.imag(gt_ker_2D)

        Din = 2*kRO*kNE1*srcCHA
        Dout = 2*dstCHA*oE1

        params = dict()

        params['verbose'] = True
        params['grappa_weights_r'] = torch.from_numpy(np.transpose(gt_ker_2D_r, (1,0))).to(torch.float32)
        params['grappa_weights_i'] = torch.from_numpy(np.transpose(gt_ker_2D_i, (1,0))).to(torch.float32)

        model = GrappaAI(Din, Dout, params)
        print(model)
        print("----------------------------------------", file=sys.stderr)

        sys.stderr.flush()

    except Exception as e:
        print("Error happened in create_grappa_ai_model", file=sys.stderr)
        print(e)

    return model

# -------------------------------------------------------------------------
def apply_grappa_ai_model(gt_dataA, model, device = torch.device('cpu')):
    '''Apply the grappa ai model
    gt_dataA: [M Din] array for kspace to recon
    '''

    recon_kspace = None

    try:
        kspace = complex_2_real(gt_dataA)

        ynew = torch.from_numpy(kspace)
        ynew = ynew.to(torch.float32)
        ynew = ynew.to(device=device) 

        model = model.to(device=device)

        model.eval()
        t0 = time.time()
        with torch.no_grad():
            recon = model(ynew)
        t1 = time.time()

        recon = recon.detach().cpu().numpy().astype(np.float32)
        recon_kspace = real_2_complex(recon)
    except Exception as e:
        print("Error happened in apply_grappa_ai_model", file=sys.stderr)
        print(e)

    return recon_kspace

# -------------------------------------------------------------------------
if __name__ == "__main__":

    GT_HOME = os.environ['GADGETRON_HOME']
    model_dir = os.path.join(GT_HOME, 'share/gadgetron/python')
    print("GT_HOME is", GT_HOME)

    GT_CMR_ML_UT_HOME = os.environ['GT_CMR_ML_UNITTEST_DIRECTORY']
    print("GT_CMR_ML_UT_HOME is", GT_CMR_ML_UT_HOME)

    gt_py_home = os.path.join(GT_HOME, 'share/gadgetron/python')
    sys.path.insert(0, gt_py_home)

    import gadgetron_toolbox_mri_core_python as gt

    data_name = os.path.join(GT_CMR_ML_UT_HOME, 'data/AI_recon/meas_MID00348_FID09272_SNR_R2/Grappa_test_data.mat')
    print(data_name)

    D = loadmat(data_name)

    ref_src = D['ref_src']
    ref_dst = D['ref_dst']
    data = D['data']

    print('ref_src = ', ref_src.shape, np.linalg.norm(ref_src))
    print('ref_dst = ', ref_dst.shape, np.linalg.norm(ref_dst))
    print('data = ', data.shape, np.linalg.norm(data))

    RO, E1, srcCHA = data.shape

    accel = 2
    oE1 = 2

    # ------------------------------------------------------

    grappa = gt.grappa2D()

    kRO = 5
    kNE1 = 4
    fitItself = True
    thres = 1e-4

    res = grappa.initialize(accel, kRO, kNE1, fitItself, thres)

    if(fitItself):
        grappa.calib(ref_src, ref_dst)
    else:
        grappa.calib(ref_src, ref_src)

    try:
        data2D = data[:,:,:,0]
        fullkspace = grappa.recon(data[:,:,:,0], True)
    except:
        fullkspace = grappa.recon(data, True)
        data2D = data
    
    if(not fitItself):
        fullkspace = fullkspace + data2D

    gt_A = grappa.get_A()
    gt_B = grappa.get_B()
    gt_ker = grappa.get_ker()

    gt_dataA = grappa.get_data_A()
    gt_dataAInd = grappa.get_data_A_index()

    try:
        kRO, kNE1, srcCHA, dstCHA, oE1 = gt_ker.shape
        gt_ker_2D = np.reshape(gt_ker, (kRO*kNE1*srcCHA, dstCHA*oE1), order='F')
    except:
        gt_ker_2D = gt_ker

    print('gt_A = ', gt_A.shape, np.linalg.norm(gt_A))
    print('gt_B = ', gt_B.shape, np.linalg.norm(gt_B))
    print('gt_ker = ', gt_ker.shape, np.linalg.norm(gt_ker))
    print('gt_dataA = ', gt_dataA.shape, np.linalg.norm(gt_dataA))
    print('gt_dataAInd = ', gt_dataAInd.shape, np.linalg.norm(gt_dataAInd))

    # ------------------------------------------------------

    if torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    print('using device:', device)

    if torch.cuda.device_count() > 1:
        print("Let's use", torch.cuda.device_count(), "GPUs!")

    model = create_grappa_ai_model(gt_ker)
    recon = apply_grappa_ai_model(gt_dataA, model, device=device)

    # get the reconed kspace
    recon_ai = grappa.recon_fill_back_kspace(recon, gt_dataAInd, RO, E1)

    output = {'grappa':fullkspace, 'ai':recon_ai}

    res_dir = os.path.join(GT_CMR_ML_UT_HOME, 'result', 'AI_recon', 'meas_MID00348_FID09272_SNR_R2')
    os.makedirs(res_dir, exist_ok=True)

    res_file = os.path.join(res_dir, 'grappa_ai_res.mat')
    print('save :', res_file)
    savemat(res_file,output)
