
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

from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils

import os
import sys
import math
import time
import random

from gadgetron_cmr_landmark_detection_util import get_landmark_from_prob, get_landmark_from_prob_fast
from gadgetron_cmr_landmark_detection_util import cpad_2d

def load_model_cmr_landmark_detection(model_dir, model_file):

    m = None

    try:
        print("Load cmr landmark model : %s" % os.path.join(model_dir, model_file), file=sys.stderr)
        t0 = time.time()
        m = torch.jit.load(os.path.join(model_dir, model_file))
        t1 = time.time()
        print("Model loading took %f seconds " % (t1-t0), file=sys.stderr)
        sys.stderr.flush()
    except Exception as e:
        print("Error happened in load_model_cmr_landmark_detection for %s" % model_file, file=sys.stderr)
        print(e)

    return m

def perform_cmr_landmark_detection_with_model(im, model_dir, model_file, p_thresh=0.1, oper_RO=352, oper_E1=352):
    model = torch.jit.load(os.path.join(model_dir, model_file))
    #print(model, file=sys.stderr)
    sys.stderr.flush()

    pts, probs = perform_cmr_landmark_detection(im, model, p_thresh=p_thresh, oper_RO=oper_RO, oper_E1=oper_E1)
    return pts, probs

def perform_cmr_landmark_detection_fast(im, model, p_thresh=0.1, oper_RO=352, oper_E1=352):
    pts, probs = perform_cmr_landmark_detection(im, model, fast_mode=1.0, p_thresh=p_thresh, oper_RO=oper_RO, oper_E1=oper_E1)

    print("perform_cmr_landmark_detection_fast, pos 1 ", file=sys.stderr)
    return pts, probs

def perform_cmr_landmark_detection(im, model, fast_mode=0.0, batch_size=8, p_thresh=0.1, oper_RO=352, oper_E1=352):
    """
    Perform CMR landmark detection

    Input :
    im : [RO, E1, N],image
    model : loaded model
    p_thres: if max(prob)<p_thres, then no landmark is found
    fast_mode : if True, use the max peak as landmark, faster
    oper_RO, oper_E1: expected array size of model. Image will be padded.

    Output:
    pts: [N_pts, 2, N], landmark points, if no landmark, it is -1
    probs: [RO, E1, N_pts+1, N], probability of detected landmark points
    """

    ori_shape = im.shape
    RO = im.shape[0]
    E1 = im.shape[1]

    try:
        im = np.reshape(im, (RO, E1, np.prod(ori_shape[2:])))
        N = im.shape[2]
    except:
        im = np.expand_dims(im, axis=2)
        RO, E1, N = im.shape

    print("perform_cmr_landmark_detection, input is ", im.shape, file=sys.stderr)
    print("perform_cmr_landmark_detection, im is ", np.linalg.norm(im), file=sys.stderr)
    print("perform_cmr_landmark_detection, fast_mode is ", fast_mode, file=sys.stderr)
    sys.stderr.flush()

    im_used, s_ro, s_e1 = cpad_2d(im, oper_RO, oper_E1)

    for n in range(N):
        #print("max for %d is %f" % (n, np.max(im_used[:,:,n])))
        im_used[:,:,n] = im_used[:,:,n] / np.max(im_used[:,:,n])

    print("perform_cmr_landmark_detection, im_used is ", np.linalg.norm(im_used), file=sys.stderr)
    sys.stderr.flush()

    im_used = np.transpose(im_used, (2,0,1))
    im_used = np.expand_dims(im_used, axis=1).astype(np.float32)
    print("perform_cmr_landmark_detection, images for model is ", im_used.shape, file=sys.stderr)
    sys.stderr.flush()

    t0 = time.time()
    images = torch.from_numpy(im_used).float()
    
    print("perform_cmr_landmark_detection, images to torch for model is ", torch.norm(images), file=sys.stderr)
    sys.stderr.flush()

    if torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    print('perform_cmr_landmark_detection, using device:', device, file=sys.stderr)

    model = model.to(device=device)
    images = images.to(device=device, dtype=torch.float32)

    #if(device == torch.device('cpu')):
    #    model.eval() 
    #    with torch.no_grad():
    #        scores = model(images)
    #else:

    ind = np.arange(0, N, int(batch_size))
    scores = torch.zeros((N, 4, images.shape[2], images.shape[3]))

    model.eval() 
    with torch.no_grad():
        for b in range(ind.shape[0]-1):
            print("processing ", ind[b], ind[b+1], file=sys.stderr)
            sys.stderr.flush()
            scores[ind[b]:ind[b+1],:,:,:] = model(images[ind[b]:ind[b+1],:,:,:]).to(device=torch.device('cpu'))

        scores[ind[-1]:N,:,:,:] = model(images[ind[-1]:N,:,:,:]).to(device=torch.device('cpu'))

    probs = torch.softmax(scores, dim=1)
    t1 = time.time()
    probs = probs.to(device=torch.device('cpu'))
    print("perform_cmr_landmark_detection, model runs in %.2f seconds " % (t1-t0))
    print("perform_cmr_landmark_detection, probs is ", probs.shape, file=sys.stderr)
    sys.stderr.flush()

    probs = probs.numpy()
    probs = probs.astype(np.float32)
    probs = np.transpose(probs, (2,3,1,0))

    print("perform_cmr_landmark_detection, probs in numpy is ", probs.shape, file=sys.stderr)
    sys.stderr.flush()

    C = probs.shape[2]

    if(s_ro>=0):
        probs = probs[s_ro:s_ro+RO, :, :, :]
    else:
        probs_used = np.reshape(probs, (probs.shape[0], probs.shape[1], C*N))
        probs_used, s_ro_p, s_e1_p = cpad_2d(probs_used, RO, probs.shape[1])
        probs = np.reshape(probs_used, (probs_used.shape[0], probs_used.shape[1], C, N))

    if(s_e1>=0):
        probs = probs[:, s_e1:s_e1+E1, :, :]
    else:
        probs_used = np.reshape(probs, (probs.shape[0], probs.shape[1], C*N))
        probs_used, s_ro_p, s_e1_p = cpad_2d(probs_used, probs.shape[0], E1)
        probs = np.reshape(probs_used, (probs_used.shape[0], probs_used.shape[1], C, N))

    print("perform_cmr_landmark_detection, probs after cpad is ", probs.shape, file=sys.stderr)
    sys.stderr.flush()

    N_pts = C-1

    pts = np.zeros((N_pts,2,N))-1.0
    print("N_pt is ", N_pts, file=sys.stderr)
    sys.stderr.flush()

    for n in range(N):
        print("perform_cmr_landmark_detection, compute landmark for image ", n, file=sys.stderr)
        sys.stderr.flush()

        for p in range(N_pts):
            prob = probs[:,:,p+1,n]
            if(fast_mode==1.0):
                pt = get_landmark_from_prob_fast(prob, thres=p_thresh)
            else:
                pt = get_landmark_from_prob(prob, thres=p_thresh, mode="mean", binary_mask=False)
            if(pt is not None):
                pts[p, 0, n] = pt[0]
                pts[p, 1, n] = pt[1]

    probs = probs.astype(np.float32)
    pts = pts.astype(np.float32)

    probs = np.reshape(probs, (RO, E1, C) + ori_shape[2:])
    pts = np.reshape(pts, (N_pts, 2) + ori_shape[2:])

    print("perform_cmr_landmark_detection, pts at final is ", pts, file=sys.stderr)
    sys.stderr.flush()
    print("perform_cmr_landmark_detection, probs at final is ", np.linalg.norm(probs), file=sys.stderr)
    sys.stderr.flush()
    print("perform_cmr_landmark_detection, probs at final is ", probs.shape, file=sys.stderr)
    sys.stderr.flush()

    print("perform_cmr_landmark_detection, pos 1 ", file=sys.stderr)
    sys.stderr.flush()
    return pts, probs

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='CMR landmark detection')
    parser.add_argument("--input", default="image.npy")
    parser.add_argument("--output", default="pts.npy")
    parser.add_argument("--prob", default="prob.npy")
    parser.add_argument("--RO", type=int, default=352)
    parser.add_argument("--E1", type=int, default=352)
    parser.add_argument("--thres", type=float, default=0.1)
    parser.add_argument("--batch_size", type=int, default=8)

    parser.add_argument("--model", default="CMR_landmark_network_ch2_ch3_ch4_myo_pts_LossMultiSoftProb_KLD_Dice_CMR_View__Pytorch_1.5.0_2020-06-06_20200606_034214.pts")

    args = parser.parse_args()

    GT_HOME = os.environ['GADGETRON_HOME']
    model_dir = os.path.join(GT_HOME, 'share/gadgetron/python')
    print("GT_HOME is", GT_HOME)

    GT_CMR_ML_UT_HOME = os.environ['GT_CMR_ML_UNITTEST_DIRECTORY']
    print("GT_CMR_ML_UT_HOME is", GT_CMR_ML_UT_HOME)

    ut_mode = True
    if(os.path.isfile(args.input)):
        ut_mode = False
        print("load - ", args.input)
        im = np.load(args.input)
        print(im.shape)

        model = load_model_cmr_landmark_detection(model_dir, args.model)
        #print(model)

        pts, probs = perform_cmr_landmark_detection(im, model, batch_size=args.batch_size, p_thresh=args.thres, oper_RO=args.RO, oper_E1=args.E1)
        print('PTs ', pts.shape)
        print('Probs ', probs.shape)

        print("save - ", args.output)
        np.save(args.output, pts)

        print("save - ", args.prob)
        np.save(args.prob, probs)

    else:

        ch2_ch3_ch4_model_file = 'CMR_landmark_network_ch2_ch3_ch4_myo_pts_LossMultiSoftProb_KLD_Dice_CMR_View__Pytorch_1.5.0_2020-06-06_20200606_034214.pts'
        sax_model_file = 'CMR_landmark_network_RO_352_E1_352_sax_with_T1_T1T2_LGE_LossMultiSoftProb_KLD_Dice_CMR_View__Pytorch_1.5.0_2020-06-15_20200615_125215.pts'

        # Unit Test
        print("=================================================================")
        print("Test RetroCine, CH4")

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 'cmr_landmark_detection', 'RetroCine', 'CH4', '20180104_1462193_ch4_ED.npy')
        print("load - ", data_file)
        im = np.load(data_file)
        print(im.shape)

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 'cmr_landmark_detection', 'RetroCine', 'CH4', '20180104_1462193_ch4_ES.npy')
        print("load - ", data_file)
        im2 = np.load(data_file)
        print(im2.shape)

        model = load_model_cmr_landmark_detection(model_dir, ch2_ch3_ch4_model_file)
        print(model)

        pts, probs = perform_cmr_landmark_detection(np.dstack((im, im2)), model, p_thresh=0.1)
        print('PTs ', pts.shape)
        print('Probs ', probs.shape)

        res_dir = os.path.join(GT_CMR_ML_UT_HOME, 'result', 'cmr_landmark_detection')
        if(os.path.isdir(res_dir)==False):
            os.mkdir(res_dir)

        res_file = os.path.join(res_dir, 'RetronCine_CH4_images.npy')
        print("save - ", res_file)
        np.save(res_file, np.dstack((im, im2)))

        res_file = os.path.join(res_dir, 'RetronCine_CH4_pts.npy')
        print("save - ", res_file)
        np.save(res_file, pts)

        res_file = os.path.join(res_dir, 'RetronCine_CH4_probs.npy')
        print("save - ", res_file)
        np.save(res_file, probs)

        # import matplotlib
        # import matplotlib.pyplot as plt
        # matplotlib.use('TkAgg')

        # plt.figure(figsize=(16,8))
        # plt.subplot(1,2,1)
        # plt.imshow(im)
        # plt.plot(pts[0,:,0], pts[0, :,1], 'r+', markersize=16)
        # plt.axis('off')
        # plt.subplot(1,2,2)
        # plt.imshow(im)
        # plt.plot(pts[1,:,0], pts[1, :,1], 'r+', markersize=16)
        # plt.axis('off')

        print("=================================================================")
        print("Test RetroCine, CH2")

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 'cmr_landmark_detection', 'RetroCine', 'CH2', '20180323_425511678_ch2_ED.npy')
        print("load - ", data_file)
        im = np.load(data_file)
        print(im.shape)

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 'cmr_landmark_detection', 'RetroCine', 'CH2', '20180323_425511678_ch2_ES.npy')
        print("load - ", data_file)
        im2 = np.load(data_file)
        print(im2.shape)

        pts, probs = perform_cmr_landmark_detection(np.dstack((im, im2)), model, p_thresh=0.1)
        print('PTs ', pts.shape)
        print('Probs ', probs.shape)

        res_dir = os.path.join(GT_CMR_ML_UT_HOME, 'result', 'cmr_landmark_detection')
        if(os.path.isdir(res_dir)==False):
            os.mkdir(res_dir)

        res_file = os.path.join(res_dir, 'RetronCine_CH2_images.npy')
        print("save - ", res_file)
        np.save(res_file, np.dstack((im, im2)))

        res_file = os.path.join(res_dir, 'RetronCine_CH2_pts.npy')
        print("save - ", res_file)
        np.save(res_file, pts)

        res_file = os.path.join(res_dir, 'RetronCine_CH2_probs.npy')
        print("save - ", res_file)
        np.save(res_file, probs)

        # plt.figure(figsize=(16,8))
        # plt.subplot(1,2,1)
        # plt.imshow(im)
        # plt.plot(pts[0,:,0], pts[0, :,1], 'r+', markersize=16)
        # plt.axis('off')
        # plt.subplot(1,2,2)
        # plt.imshow(im)
        # plt.plot(pts[1,:,0], pts[1, :,1], 'r+', markersize=16)
        # plt.axis('off')

        print("=================================================================")
        print("Test RTCine, CH2")

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 'cmr_landmark_detection', 'RTCine', 'CH2', 'lax_images_for_AI.npy')
        print("load - ", data_file)
        im = np.load(data_file)
        print(im.shape)

        pts, probs = perform_cmr_landmark_detection(im, model, p_thresh=0.1, fast_mode=1)
        print('PTs ', pts.shape)
        print('Probs ', probs.shape)

        res_dir = os.path.join(GT_CMR_ML_UT_HOME, 'result', 'cmr_landmark_detection')
        if(os.path.isdir(res_dir)==False):
            os.mkdir(res_dir)

        res_file = os.path.join(res_dir, 'RTCine_CH2_images.npy')
        print("save - ", res_file)
        np.save(res_file, im)

        res_file = os.path.join(res_dir, 'RTCine_CH2_pts.npy')
        print("save - ", res_file)
        np.save(res_file, pts)

        res_file = os.path.join(res_dir, 'RTCine_CH2_probs.npy')
        print("save - ", res_file)
        np.save(res_file, probs)

        print("=================================================================")
        print("Test RTCine, CH4")

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 'cmr_landmark_detection', 'RTCine', 'CH4', 'meas_MID00176_FID134706_2ch_RT_cine_160res_stack_Physio_Interp.npy')
        print("load - ", data_file)
        im = np.load(data_file)
        print(im.shape)

        pts, probs = perform_cmr_landmark_detection(im, model, p_thresh=0.1, fast_mode=1)
        print('PTs ', pts.shape)
        print('Probs ', probs.shape)

        res_dir = os.path.join(GT_CMR_ML_UT_HOME, 'result', 'cmr_landmark_detection')
        if(os.path.isdir(res_dir)==False):
            os.mkdir(res_dir)

        res_file = os.path.join(res_dir, 'RTCine_CH4_images.npy')
        print("save - ", res_file)
        np.save(res_file, im)

        res_file = os.path.join(res_dir, 'RTCine_CH4_pts.npy')
        print("save - ", res_file)
        np.save(res_file, pts)

        res_file = os.path.join(res_dir, 'RTCine_CH4_probs.npy')
        print("save - ", res_file)
        np.save(res_file, probs)

        print("=================================================================")
        print("Test RetroCine, T1")

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 't1t2', 'T1SR_Mapping_SASHA_HC_T1T2_42363_622647050_622647055_2914_20200506-143306', 't1.npy')
        print("load - ", data_file)
        im = np.load(data_file)
        print(im.shape)

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 't1t2', 'T1SR_Mapping_SASHA_HC_T1T2_42363_622647050_622647055_2914_20200506-143306', 't2.npy')
        print("load - ", data_file)
        im2 = np.load(data_file)
        print(im2.shape)

        model = load_model_cmr_landmark_detection(model_dir, sax_model_file)
        print(model)

        pts, probs = perform_cmr_landmark_detection(np.dstack((im, im2)), model, p_thresh=0.1)
        print('PTs ', pts.shape)
        print('Probs ', probs.shape)

        res_file = os.path.join(res_dir, 't1t2_SAX_images.npy')
        print("save - ", res_file)
        np.save(res_file, np.dstack((im, im2)))

        res_file = os.path.join(res_dir, 't1t2_SAX_pts.npy')
        print("save - ", res_file)
        np.save(res_file, pts)

        res_file = os.path.join(res_dir, 't1t2_SAX_probs.npy')
        print("save - ", res_file)
        np.save(res_file, probs)

        print("=================================================================")
        print("Test RetroCine, T1T2, T1, inline")

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 't1t2', 'T1SR_Mapping_SASHA_HC_T1T2_42363_665654346_665654351_47_20200609-100835', 't1.npy')
        print("load - ", data_file)
        im = np.load(data_file)
        print(im.shape)

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 't1t2', 'T1SR_Mapping_SASHA_HC_T1T2_42363_665654346_665654351_47_20200609-100835', 't2.npy')
        print("load - ", data_file)
        im2 = np.load(data_file)
        print(im2.shape)

        pts, probs = perform_cmr_landmark_detection_with_model(im2, model_dir, sax_model_file, p_thresh=0.1)

        #model = load_model_cmr_landmark_detection(model_dir, 'CMR_landmark_network_sax_with_T1_LossMultiSoftProb_KLD_Dice_CMR_View__Pytorch_1.5.0_2020-06-06_20200606_204416.pts')
        #print(model)

        #pts, probs = perform_cmr_landmark_detection(im2, model, p_thresh=0.1)
        print('PTs ', pts.shape)
        print('Probs ', probs.shape)

        res_dir = os.path.join(GT_CMR_ML_UT_HOME, 'result', 'cmr_landmark_detection', 'T1SR_Mapping_SASHA_HC_T1T2_42363_665654346_665654351_47_20200609-100835')
        if(os.path.isdir(res_dir)==False):
            os.mkdir(res_dir)

        res_file = os.path.join(res_dir, 't1t2_SAX_images.npy')
        print("save - ", res_file)
        np.save(res_file, im2)

        res_file = os.path.join(res_dir, 't1t2_SAX_pts.npy')
        print("save - ", res_file)
        np.save(res_file, pts)

        res_file = os.path.join(res_dir, 't1t2_SAX_probs.npy')
        print("save - ", res_file)
        np.save(res_file, probs)

        print("=================================================================")
        print("Test RetroCine, T1T2, T1, inline 2")

        data_file = os.path.join(GT_CMR_ML_UT_HOME, 'data', 't1t2', 'T1SR_Mapping_SASHA_HC_T1T2_42363_665654346_665654351_47_20200609-100835', 'map_t2_encoding_0_processing_4_resized.npy')
        print("load - ", data_file)
        im = np.load(data_file)
        print(im.shape)

        im = np.transpose(im, (1, 0))

        pts, probs = perform_cmr_landmark_detection_with_model(im, model_dir, sax_model_file, p_thresh=0.1)

        #model = load_model_cmr_landmark_detection(model_dir, 'CMR_landmark_network_sax_with_T1_LossMultiSoftProb_KLD_Dice_CMR_View__Pytorch_1.5.0_2020-06-06_20200606_204416.pts')
        #print(model)

        #pts, probs = perform_cmr_landmark_detection(im2, model, p_thresh=0.1)
        print('PTs ', pts.shape)
        print('Probs ', probs.shape)

        res_dir = os.path.join(GT_CMR_ML_UT_HOME, 'result', 'cmr_landmark_detection', 'T1SR_Mapping_SASHA_HC_T1T2_42363_665654346_665654351_47_20200609-100835')
        if(os.path.isdir(res_dir)==False):
            os.mkdir(res_dir)

        res_file = os.path.join(res_dir, 't1t2_SAX2_images.npy')
        print("save - ", res_file)
        np.save(res_file, im2)

        res_file = os.path.join(res_dir, 't1t2_SAX2_pts.npy')
        print("save - ", res_file)
        np.save(res_file, pts)

        res_file = os.path.join(res_dir, 't1t2_SAX2_probs.npy')
        print("save - ", res_file)
        np.save(res_file, probs)
