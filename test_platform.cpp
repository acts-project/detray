// test_macros.cpp
#if defined(__HIP_PLATFORM_AMD__)
#pragma message("✅ [HIP] AMD platform is active.")
#elif defined(__HIP_PLATFORM_NVCC__)
#pragma message("✅ [HIP] NVCC (NVIDIA) platform is active.")
#else
#pragma message("❌ No HIP platform macro defined.")
#endif

#if defined(__HIPCC__)
#pragma message("✅ [HIP] platform is active.")
#else
#pragma message("✅ [HIP] platform is active.")
#endif

int main() {
    return 0;
}
