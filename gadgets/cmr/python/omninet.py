#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OmniNet inference wrapper for Gadgetron's CmrOmniGadget.
"""

import os
import gc
import time

import numpy as np
import torch

# Scan types this deployment supports. Static LGE/DB/WB are refused: they cannot be told apart from
# the header's encoding structure (all three have n_phase = n_rep = n_segment = 1). The artifact
# itself can reconstruct them -- only the header-based identification is impossible.
DEPLOY_FILE_TYPES = ("Perfusion", "Retro", "RT")

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

# ---------------------------------------------------------

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
        torch.set_float32_matmul_precision("high")
        torch._C._jit_override_can_fuse_on_gpu(False)
        torch._C._jit_set_texpr_fuser_enabled(False)
        model = torch.jit.load(model_file_name, map_location=map_location).eval()
        logger.info(f"---> Model loading took {time.time()-t0:.2f}s")
    except Exception as e:
        logger.error(f"Error in load_model_for_inference for {model_file_name}")
        logger.exception(e)
        model = None
    return model

# ---------------------------------------------------------

def infer(model, kspace: torch.Tensor, sens_maps: torch.Tensor, file_type: str, NumOfProtonDensityImages: None | int) -> torch.Tensor:
    """Full-depth reconstruction.

    kspace:    (slice, coils, depth, ro, pe) complex64
    sens_maps: (slice, coils, ro, pe)        complex64
    file_type: str -- scan type, must be one of DEPLOY_FILE_TYPES, "Perfusion", "Retro", "RT"
    returns:   (slice, 1, depth, ro, pe)     complex64
    """
    assert kspace.ndim == 5, f"expected 5D kspace, got {tuple(kspace.shape)}"
    assert sens_maps.ndim == 4, f"expected 4D sens_maps, got {tuple(sens_maps.shape)}"

    if file_type not in DEPLOY_FILE_TYPES:
        raise ValueError(f"unsupported scan ({file_type or 'static LGE/DB/WB'})")

    # ft_idx comes from the ARTIFACT's own ordering, never hardcoded: it is a bare int, so a stale
    # hardcoded order would select the wrong scan type SILENTLY.
    ft_idx = list(model.file_types).index(file_type)

    is_pd = None
    if file_type == "Perfusion" and NumOfProtonDensityImages is not None:
        D = kspace.shape[2]
        is_pd = torch.arange(D, device=kspace.device) < (NumOfProtonDensityImages or 0)

    with torch.inference_mode():
        # the scripted model derives the sampling mask from the k-space support itself
        return model(kspace, sens_maps, ft_idx, is_pd)

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
    if arr_roe1cha.ndim == 5:
        # (RO, E1, CHA, N, SLC) -> torch (SLC, CHA, N, RO, E1) on `device`.
        x = np.ascontiguousarray(arr_roe1cha.transpose(4, 2, 3, 0, 1))
    else:
        # (RO, E1, CHA, SLC) -> torch (SLC, CHA, RO, E1) on `device`.
        x = np.ascontiguousarray(arr_roe1cha.transpose(3, 2, 0, 1))
    return torch.from_numpy(x).to(dtype=dtype, device=device)

# ---------------------------------------------------------

def apply_omninet(kspace, coil_map, model, file_type="RT", NumOfProtonDensityImages=None,
                  use_cpu=False,
                  verbose=True):
    """Run OmniNet once per batch element and assemble the output.

    kspace: `(ro, pe, coils, depth, slice)` cfloat k-space
    sens_maps: `(ro, pe, coils, slice)` cfloat sensitivity maps
    """
    RO, E1, CHA, N, SLC = kspace.shape
    if coil_map.shape != (RO, E1, CHA, SLC):
        raise ValueError(f"apply_omninet: coil_map shape {coil_map.shape} != expected {(RO, E1, CHA, SLC)}")

    if verbose:
        logger.info(f"---> apply_omninet kspace  {kspace.shape} {kspace.dtype}")
        logger.info(f"---> apply_omninet coil_map  {coil_map.shape} {coil_map.dtype}")
        logger.info(f"---> apply_omninet file_type  {file_type}")
        logger.info(f"---> apply_omninet NumOfProtonDensityImages  {NumOfProtonDensityImages}")
        logger.info(f"---> apply_omninet RO={RO} E1={E1} CHA={CHA} N={N} SLC={SLC}")

    # Always allocate the output up front so the error path returns a valid shape.
    output = np.zeros((RO, E1, CHA, N, SLC), dtype=np.csingle, order="F")

    try:
        t0 = time.time()
        if use_cpu:
            device = torch.device("cpu")
            model.eval()
        else:
            model, devices = prep_model(model)
            device = devices[0]

        with torch.inference_mode():
            kspace_tensor = _to_model_input(kspace, device, torch.cfloat)
            coil_map_tensor = _to_model_input(coil_map, device, torch.cfloat)
            out = infer(model, kspace_tensor, coil_map_tensor, file_type, NumOfProtonDensityImages)
            output = out.detach().to(dtype=torch.cfloat, device="cpu").numpy()
            output = np.asfortranarray(output.transpose(3, 4, 1, 2, 0))  # (SLC, 1, N, RO, E1) -> (RO, E1, 1, N, SLC)

        logger.info(f"---> apply_omninet took {time.time()-t0:.3f}s for kspace {kspace.shape}")
    except Exception as e:
        logger.error("Error in apply_omninet; returning zeros")
        logger.exception(e)
        output = np.zeros((RO, E1, CHA, N, SLC), dtype=np.csingle, order="F")
    finally:
        finish_model()

    logger.info(f"---> apply_omninet output {output.shape} {output.dtype}")
    return output

# ---------------------------------------------------------
if __name__ == "__main__":
    pass