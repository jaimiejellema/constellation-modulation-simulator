//
// Created by jaimie on 22-5-19.
//

#ifndef GAM_CONSTANTS_H
#define GAM_CONSTANTS_H

#ifdef __JETBRAINS_IDE__
    #define __host__
    #define __device__
    #define __shared__
    #define __constant__
    #define __global__

    // This is slightly mental, but gets it to properly index device function calls like __popc and whatever.
    #define __CUDACC__
    #include <device_functions.h>

    // These headers are all implicitly present when you compile CUDA with clang. Clion doesn't know that, so
    // we include them explicitly to make the indexer happy. Doing this when you actually build is, obviously,
    // a terrible idea :D
    #include <__clang_cuda_builtin_vars.h>
    #include <__clang_cuda_intrinsics.h>
    #include <__clang_cuda_math_forward_declares.h>
    #include <__clang_cuda_complex_builtins.h>
    #include <__clang_cuda_cmath.h>
#endif // __JETBRAINS_IDE__

#endif //GAM_CONSTANTS_H
