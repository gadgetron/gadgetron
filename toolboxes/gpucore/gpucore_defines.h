#pragma once

#if defined(__GNUC__)
        
#define __gannotate__(a) \
        __attribute__((a))

#define __glocation__(a) \
        __gannotate__(a)

#elif defined(_WIN32)

#define __inline__ \
        __inline

#define __gannotate__(a) \
        __declspec(a)

#define __glocation__(a) \
        __gannotate__(__##a##__)

#endif

#if !defined(__CUDACC__) && !defined(__CUDABE__)

#undef __gannotate__
#define __gannotate__(a)

#endif

#define __gad_host__ \
        __glocation__(host)

#define __gad_device__ \
        __glocation__(device)
