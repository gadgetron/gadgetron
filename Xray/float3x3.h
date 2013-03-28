#pragma once

//#define SINCOSF __sincosf // fast math
#define SINCOSF sincosf
#include "vector_td.h"

namespace Gadgetron{

struct float3x3 {
    floatd3 row0;
    floatd3 row1;
    floatd3 row2;
};

__inline__ __host__ __device__ 
float3x3 make_float3x3(float v0, float v1, float v2,
                       float v3, float v4, float v5,
                       float v6, float v7, float v8) {
    float3x3 m;
    m.row0 = floatd3(v0, v1, v2);
    m.row1 = floatd3(v3, v4, v5);
    m.row2 = floatd3(v6, v7, v8);
    return m;
}

__inline__ __device__ 
floatd3 mul(float3x3 m, floatd3 v) {
    return floatd3( dot(m.row0,v), dot(m.row1,v), dot(m.row2,v) );
}


__inline__ __device__ float3x3 calcRotationMatrixAroundX(float angle) {
    float cos_angle, sin_angle;
    SINCOSF(angle, &sin_angle, &cos_angle);
  
    // Build projection rotation matrix
    float3x3 rotation = make_float3x3(1,         0,          0,
                                      0, cos_angle, -sin_angle,
                                      0, sin_angle,  cos_angle);
    return rotation;
}

__inline__ __device__ float3x3 calcRotationMatrixAroundY(float angle) {
    float cos_angle, sin_angle;
    SINCOSF(angle, &sin_angle, &cos_angle);
  
    // Build projection rotation matrix
    float3x3 rotation = make_float3x3( cos_angle, 0, sin_angle,
                                               0, 1,         0,
                                      -sin_angle, 0, cos_angle);
    return rotation;
}

__inline__ __device__ float3x3 calcRotationMatrixAroundZ(float angle) {
    float cos_angle, sin_angle;
    sincosf(angle, &sin_angle, &cos_angle);
  
    // Build projection rotation matrix
    float3x3 rotation = make_float3x3(cos_angle, -sin_angle, 0,
                                      sin_angle,  cos_angle, 0,
                                              0,          0, 1);
    return rotation;
}


}
