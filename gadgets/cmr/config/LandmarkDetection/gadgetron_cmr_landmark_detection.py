
"""Perform cmr landmark detection"""

import os
import time

import numpy as np
from scipy import interpolate
from scipy.special import softmax

import onnxruntime as ort

from gadgetron_cmr_landmark_detection_util import get_landmark_from_prob, cpad_2d

def load_model_onnx(model_dir, model_file):
    """Load onnx format model

    Args:
        model_dir (str): folder to store the model; if None, only model file is used
        model_file (str): model file name, can be the full path to the model
        
    Returns:
        model: loaded model
    """
    
    m = None
    
    try:
        if(model_dir is not None):
            model_full_file = os.path.join(model_dir, model_file)
        else:
            model_full_file = model_file

        print("gadgetron_cmr_landmark_detection, load model : %s" % model_full_file)
        t0 = time.time()       
        m = ort.InferenceSession(model_full_file, providers=['CPUExecutionProvider'])
        t1 = time.time()
        print(f"gadgetron_cmr_landmark_detection, model loading took {t1-t0:.2f} seconds ")
    except Exception as e:
        print("Error happened in load_model_onnx for %s" % model_file)
        raise e

    return m

def find_landmark_from_prob_maps(probs, p_thres=0.1, mode='mean'):
    """
    find landmark points from probs

    Input:
        probs : [RO, E1, C, N]

    Output:

        pts: [C, 2, N], landmar points, -1 means no landmark
    """

    RO, E1, C, N = probs.shape
    N_pts = C

    pts = np.zeros((N_pts,2,N))-1.0
    print(f"N_pt is {N_pts}")

    for n in range(N):
        print(f"find_landmark_from_prob_maps, compute landmark for image {n}")

        for p in range(N_pts):
            prob = probs[:,:,p,n]
            pt = get_landmark_from_prob(prob, thres=p_thres, mode="mean", binary_mask=False)
            if(pt is not None):
                pts[p, 0, n] = pt[0]
                pts[p, 1, n] = pt[1]

    pts = pts.astype(np.float32)
    print(f"find_landmark_from_prob_maps, pts at final is {pts}")

    return pts


def fill_missing_pts(pts, smoothing=1, s=-1):
    """
    Fill in missing landmarks using interpolation

    Input: 
        pts: [Npts, D, PHS, SLC], missing landmarks are -1
        smoothing: if 1, output smoothed pts; if 0, only replace missing pts
        s : the amount of smoothing to perform during the spline fit, s=0, no smoothing; if s==-1, use default smoothing (m-sqrt(2*m))
    Output:
        pts_filled: [Npts, D, PHS, SLC]
    """

    pts_filled = pts
    pts_pts_smoothed = pts

    try:
        if(pts.ndim<3):
            print(f"fill_missing_pts, wrong input pts shape {pts.shape}")
            return pts

        if(pts.ndim==3):
            pts = np.expand_dims(pts, axis=3)
            Npts, D, PHS, SLC = pts.shape

        if(pts.ndim==4):
            Npts, D, PHS, SLC = pts.shape

        if(s==-1):
            s = PHS - np.sqrt(2*PHS)

        pts_filled = pts.copy()
        pts_smoothed = pts.copy()

        for slc in range(SLC):
            for p in range(Npts):
                x = pts[p, 0, :, slc].squeeze().copy()
                y = pts[p, 1, :, slc].squeeze().copy()

                ind_x = np.argwhere(x > 0)
                x_used = x[x > 0].astype(np.float64) 
                y_used = y[x > 0].astype(np.float64) 

                x_tck = interpolate.splrep(ind_x, x_used, k=5, s=int(s))
                y_tck = interpolate.splrep(ind_x, y_used, k=5, s=int(s))
                xnew = interpolate.splev(np.arange(0, PHS), x_tck, der=0)
                ynew = interpolate.splev(np.arange(0, PHS), y_tck, der=0)

                pts_smoothed[p, 0, :, slc] = xnew
                pts_smoothed[p, 1, :, slc] = ynew

                if(smoothing):
                    pts_filled[p, 0, :, slc] = xnew
                    pts_filled[p, 1, :, slc] = ynew
                else:
                    missing_ind = np.argwhere(x <= 0)
                    if(missing_ind.shape[0]>0):
                        x[missing_ind] = xnew[missing_ind]
                        pts_filled[p, 0, :, slc] = x
                        y[missing_ind] = ynew[missing_ind]
                        pts_filled[p, 1, :, slc] = y
    except Exception as e:
        print("Error happened in fill_missing_pts ... ")
        pts_filled = pts
        pts_pts_smoothed = pts
        print(e)

    return pts_filled.astype(np.float32), pts_smoothed.astype(np.float32)

def prepare_image_for_inference(im_used):
    """Prepare input images for inference.

    Args:
        im_used ([RO, E1, N]): input images

    Returns:
        im_ready ([N, 1, RO, E1]) : images ready for inference
    """

    RO, E1, N = im_used.shape

    t0 = time.time()
    for n in range(N):
        im_used[:,:,n] = im_used[:,:,n] / np.max(im_used[:,:,n])

    print(f"cmr_landmark_detection_onnx, prepare_image_for_inference,  im_used is {np.linalg.norm(im_used)}")

    im_ready = np.transpose(im_used, (2,0,1))
    im_ready = np.expand_dims(im_ready, axis=1).astype(np.float32)
    print(f"cmr_landmark_detection_onnx, prepare_image_for_inference,  images for model is {im_ready.shape}")

    t1 = time.time()
    print(f"cmr_landmark_detection_onnx, prepare_image_for_inference,  prepare input data in {t1-t0} seconds ")

    return im_ready


def call_model_inference(im_used, model):
    """call model on input images

    Args:
        im_used ([B, 1, oper_RO, oper_E1]): images ready for processing
        model : pre-loaded model

    Returns:
        probs, [oper_RO, oper_E1, C, B]: detected probabilties
    """

    print(f"call_model_inference, perform inference on onnx model ... ")

    t0 = time.time()
    input_name = model.get_inputs()[0].name
    output_name = model.get_outputs()[0].name
    y_pred = model.run([output_name], {input_name: im_used})[0]
    probs = softmax(y_pred, axis=1)
    t1 = time.time()
    print(f"Onnx model runs in {t1 - t0} seconds ")

    probs = probs.astype(np.float32)
    probs = np.transpose(probs, (2,3,1,0))

    return probs


def get_pts_from_probs(probs, p_thresh=0.75):
    """Extract points from probablity maps

    Args:
        probs ([RO, E1, C, N]): probablities for C-1 points
        p_thresh (float, optional): threshold to extract points. Defaults to 0.75.
    """
    RO, E1, C, N = probs.shape
    N_pts = C-1

    pts = np.zeros((N_pts,2,N))-1.0
    print(f"N_pt is {N_pts}")

    t0 = time.time()
    for n in range(N):
        print(f"get_pts_from_probs, compute landmark for image {n}")

        for p in range(N_pts):
            prob = probs[:,:,p+1,n]
            pt = get_landmark_from_prob(prob, thres=p_thresh, mode="mean", binary_mask=False)
            if(pt is not None):
                pts[p, 0, n] = pt[0]
                pts[p, 1, n] = pt[1]

    t1 = time.time()
    print(f"get_pts_from_probs, compute landmark from prob in {t1-t0} seconds ")

    return pts.astype(np.float32)

def perform_cmr_landmark_detection(im, model, p_thresh=0.1, oper_RO=320, oper_E1=320):
    """
    Perform CMR landmark detection

    Input :
        im : [RO, E1, N],image
        model : loaded model
        p_thresh: if max(prob)<p_thresh, then no landmark is found
        oper_RO, oper_E1: expected array size of model. Image will be padded.

    Output:
        pts: [N_pts, 2, N], landmark points, if no landmark, it is -1
        probs: [RO, E1, N_pts+1, N], probability of detected landmark points
        im_used: [oper_RO, oper_E1, N], images used for detection
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

    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  input is {im.shape}")
    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  p_thresh is {p_thresh}")
    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  oper_RO is {oper_RO}")
    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  oper_E1 is {oper_E1}")

    im_used, _, _ = cpad_2d(im.copy(), oper_RO, oper_E1)

    im_used = prepare_image_for_inference(im_used)
    probs = call_model_inference(im_used, model)
    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  probs in numpy is {probs.shape}")

    C = probs.shape[2]
    N_pts = C-1

    probs_used = np.reshape(probs, (probs.shape[0], probs.shape[1], C*N))
    probs_used, _, _ = cpad_2d(probs_used, RO, E1)
    probs = np.reshape(probs_used, (probs_used.shape[0], probs_used.shape[1], C, N))

    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  probs after cpad is {probs.shape}")

    pts = get_pts_from_probs(probs, p_thresh=p_thresh)

    probs = np.reshape(probs, (RO, E1, C) + ori_shape[2:])
    pts = np.reshape(pts, (N_pts, 2) + ori_shape[2:])

    if(pts.ndim==2):
        pts = np.expand_dims(pts, axis=2)

    im_used = np.squeeze(im_used)
    if(im_used.ndim==3):
        im_used = np.transpose(im_used, [1, 2, 0])

    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  pts at final is {pts}")
    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  probs at final is {np.linalg.norm(probs)}")
    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  probs at final is {probs.shape}")
    print(f"gadgetron_cmr_landmark_detection, perform_cmr_landmark_detection,  im_used at final is {im_used.shape}")

    return pts.astype(np.float32), probs.astype(np.float32), im_used.astype(np.float32)