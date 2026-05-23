#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PSIRNet inference wrapper for Gadgetron's CmrPSIRNetGadget.

Shape contract with the gadget (CmrPSIRNetGadget.cpp, perform_psir):
    IR        np.ndarray  shape (RO, E1, CHA, B) complex64
    PD        np.ndarray  shape (RO, E1, CHA, B) complex64
    coil_map  np.ndarray  shape (RO, E1, CHA, B) complex64
    model     torch.jit.ScriptModule with forward(a, b, c, d) ->
                  a, b, d : (1, CHA, RO, E1)  complex
                  c       : (1, 1,   RO, E1)  bool
              returning (1, 1, RO, E1) real32 (or anything squeezable to RO,E1).

PSIRNet is a single-shot model (one IR + one PD k-space per call), see
psirnet/src/models.py:PSIRNet docstring. The C++ caller packs every
(average, slice) pair into the flat batch dim `B = N*SLC` (as `b = n + slc*N`)
and unpacks the output the same way. We run exactly one inference per
batch element and stay agnostic to what the batch represents.

Returns:
    psir      np.ndarray  shape (RO, E1, 1, B) float32, Fortran-contiguous.
              PSIR is genuinely real-valued (the model emits a signed scalar
              field), so we return float32 and let the C++ side wrap into
              std::complex<float> on the way into IsmrmrdImageArray::data_.
"""

import os
import gc
import time

import numpy as np
import torch

# Optional Gadgetron helpers; fall back to stdlib if not on sys.path.
try:
    from gadgetron_logger import logger
except ImportError:
    import logging
    logger = logging.getLogger("psirnet")
    if not logger.handlers:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s %(message)s"))
        logger.addHandler(h)
    logger.setLevel(logging.INFO)

try:
    from colorama import Fore, Style
except ImportError:
    class _NoColor:
        def __getattr__(self, _): return ""
    Fore = Style = _NoColor()

CUDA_RAM_LIMIT = 8.0   # GB


# ---------------------------------------------------------

def detect_device(total_ram_threshold_in_GB=CUDA_RAM_LIMIT):
    if not torch.cuda.is_available():
        logger.info(f"{Fore.GREEN}detect_device: no CUDA, using CPU{Style.RESET_ALL}")
        return [torch.device("cpu")]

    t0 = time.time()
    visible = os.getenv("CUDA_VISIBLE_DEVICES")
    if visible:
        devices = [torch.device(f"cuda:{i}") for i in visible.split(",")]
    else:
        devices = [torch.device("cuda:0")]

    devices_valid = []
    for device in devices:
        try:
            with torch.cuda.device(device):
                free_mem, total_mem = torch.cuda.mem_get_info()
            total_gb = total_mem / 1024.0 ** 3
            free_gb = free_mem / 1024.0 ** 3
            logger.info(f"{Fore.YELLOW}GPU {device.index}{Style.RESET_ALL} "
                        f"total RAM {total_gb:.2f} GB, free {free_gb:.2f} GB")
            if total_gb >= total_ram_threshold_in_GB:
                devices_valid.append(device)
        except Exception as e:
            logger.info(f"Check CUDA RAM failed for GPU {device.index}: {e}")

    devices = devices_valid if devices_valid else [torch.device("cpu")]
    logger.info(f"detect_device took {time.time()-t0:.2f}s, using {devices}")
    return devices


# ---------------------------------------------------------

def load_model_for_inference(model_dir, model_file):
    """Load a TorchScript model from <model_dir>/<model_file>.

    Loads directly onto CUDA when available so that frozen weights AND
    TorchScript CONSTANTS land on the same device the gadget will use for
    inputs (avoids "Input type cuda / weight type cpu" RuntimeErrors that
    nn.Module.to() can't fix because it doesn't move ScriptModule
    CONSTANTS).
    """
    model_file_name = os.path.join(model_dir, model_file) if model_dir else model_file
    map_location = "cuda:0" if torch.cuda.is_available() else "cpu"
    try:
        logger.info(f"---> Load model {model_file_name} (map_location={map_location})")
        t0 = time.time()
        model = torch.jit.load(model_file_name, map_location=map_location)
        logger.info(f"---> Model loading took {time.time()-t0:.2f}s")
    except Exception as e:
        logger.error(f"Error in load_model_for_inference for {model_file_name}")
        logger.exception(e)
        model = None
    return model


# ---------------------------------------------------------

def prep_model(model):
    if torch.get_num_threads() < os.cpu_count() * 0.5:
        torch.set_num_threads(int(os.cpu_count() * 0.8))
        logger.info(f"---> set number of cpu threads {torch.get_num_threads()}")

    devices = detect_device(total_ram_threshold_in_GB=CUDA_RAM_LIMIT)

    if isinstance(model, (torch.nn.Module, torch.jit.ScriptModule)):
        model.eval()
        # Move model to the selected device so its weights live on the same
        # device as the inputs created later in apply_psirnet().
        try:
            model = model.to(devices[0])
        except Exception as e:
            logger.warning(f"prep_model: model.to({devices[0]}) failed: {e}")
    else:
        logger.info("model is not a torch Module / ScriptModule")
    return model, devices


# ---------------------------------------------------------

def finish_model():
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


# ---------------------------------------------------------

def _to_model_input(arr_roe1cha, device, dtype):
    """(RO, E1, CHA) -> torch (1, CHA, RO, E1) on `device`."""
    x = np.ascontiguousarray(arr_roe1cha.transpose(2, 0, 1))
    return torch.from_numpy(x).to(dtype=dtype, device=device).unsqueeze(0)


def apply_psirnet(IR, PD, coil_map, model,
                  total_ram_threshold_in_GB=CUDA_RAM_LIMIT,
                  use_cpu=False,
                  verbose=True):
    """Run PSIRNet once per batch element and assemble the output.

    The 4th dim is a flat batch; the C++ caller packs (n, slc) into it as
    b = n + slc*N so the same wrapper serves both single-shot-per-slice and
    per-average-per-slice modes without code changes here.
    """
    RO, E1, CHA, B = IR.shape
    if PD.shape != IR.shape:
        raise ValueError(f"apply_psirnet: PD shape {PD.shape} != IR shape {IR.shape}")
    if coil_map.shape != IR.shape:
        raise ValueError(
            f"apply_psirnet: coil_map shape {coil_map.shape} != IR shape {IR.shape}")

    if verbose:
        logger.info(f"---> apply_psirnet IR  {IR.shape} {IR.dtype}")
        logger.info(f"---> apply_psirnet PD  {PD.shape} {PD.dtype}")
        logger.info(f"---> apply_psirnet sm  {coil_map.shape} {coil_map.dtype}")
        logger.info(f"---> apply_psirnet RO={RO} E1={E1} CHA={CHA} B={B}")

    # Always allocate the output up front so the error path returns a valid shape.
    output = np.zeros((RO, E1, 1, B), dtype=np.float32, order="F")

    try:
        t0 = time.time()
        if use_cpu:
            device = torch.device("cpu")
            model.eval()
        else:
            model, devices = prep_model(model)
            device = devices[0]

        with torch.inference_mode():
            for b in range(B):
                ir_k = IR[:, :, :, b]            # (RO, E1, CHA) one batch element
                pd_k = PD[:, :, :, b]
                sm_k = coil_map[:, :, :, b]

                a = _to_model_input(ir_k, device, torch.cfloat)
                a_pd = _to_model_input(pd_k, device, torch.cfloat)
                d = _to_model_input(sm_k, device, torch.cfloat)

                # Sampling mask from coil 0 of IR, same logic as training:
                #   mask = (ir_kspace != 0)[0:1, ...]   (psirnet/src/data.py:74)
                # Training axis order was (CHA, RO, E1); gadget gives us (RO, E1, CHA),
                # so coil-0 selection becomes ir_k[:, :, 0].
                mask_k = (ir_k[:, :, 0] != 0)      # (RO, E1) bool
                c = torch.from_numpy(np.ascontiguousarray(mask_k)).to(
                    dtype=torch.bool, device=device).unsqueeze(0).unsqueeze(0)  # (1,1,RO,E1)

                out_k = model(a, a_pd, c, d)       # expected (1, 1, RO, E1) real32
                out_k = out_k.detach().to(dtype=torch.float32, device="cpu").numpy()
                out_k = np.squeeze(out_k)
                if out_k.shape != (RO, E1):
                    raise RuntimeError(
                        f"PSIRNet output has shape {out_k.shape}, expected ({RO},{E1})")

                output[:, :, 0, b] = out_k

        logger.info(f"---> apply_psirnet took {time.time()-t0:.3f}s for B={B} batch elements")
    except Exception as e:
        logger.error("Error in apply_psirnet; returning zeros")
        logger.exception(e)
        output = np.zeros((RO, E1, 1, B), dtype=np.float32, order="F")
    finally:
        finish_model()

    logger.info(f"---> apply_psirnet output {output.shape} {output.dtype}")
    return np.asfortranarray(output)


# ---------------------------------------------------------

if __name__ == "__main__":
    # Minimal smoke test. Adjust paths to your install.
    GT_HOME = os.environ.get("GADGETRON_HOME", "/usr/local")
    model_dir = os.path.join(GT_HOME, "share/gadgetron/python/cmr_ml/models")
    model_file = "PSIRNet.pts"

    m = load_model_for_inference(model_dir, model_file)
    if m is None:
        raise SystemExit("Failed to load model")

    # Fake one batch element: 1 batch (== 1 slice, 1 average).
    RO, E1, CHA, B = 256, 192, 30, 1
    rng = np.random.default_rng(0)
    IR = np.asfortranarray(
        (rng.standard_normal((RO, E1, CHA, B)) +
         1j * rng.standard_normal((RO, E1, CHA, B))).astype(np.complex64))
    PD = np.asfortranarray(IR.copy())
    SM = np.asfortranarray(
        (rng.standard_normal((RO, E1, CHA, B)) +
         1j * rng.standard_normal((RO, E1, CHA, B))).astype(np.complex64))

    out = apply_psirnet(IR, PD, SM, m)
    print("output shape:", out.shape, "dtype:", out.dtype,
          "F-contig:", out.flags["F_CONTIGUOUS"])
