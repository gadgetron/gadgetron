#ifndef GADGETRON_CONFIG_H
#define GADGETRON_CONFIG_H

#define GADGETRON_VERSION_MAJOR 4
#define GADGETRON_VERSION_MINOR 2
#define GADGETRON_VERSION_PATCH 1
#define GADGETRON_VERSION_STRING "4.2.1"
#define GADGETRON_CONFIG_PATH "share/gadgetron/config"
#define GADGETRON_PYTHON_PATH "share/gadgetron/python"
#define GADGETRON_GIT_SHA1_HASH "54fa177320e95b20b4d676e2602fa51b74ecacf3"
#define GADGETRON_CUDA_NVCC_FLAGS " -arch=sm_80 -gencode arch=compute_80,code=sm_80  -gencode arch=compute_80,code=sm_80  -gencode arch=compute_80,code=sm_80  -gencode arch=compute_80,code=sm_80  --std=c++17"
#define GADGETRON_VAR_DIR "/var/lib/gadgetron"

#endif //GADGETRON_CONFIG_H
