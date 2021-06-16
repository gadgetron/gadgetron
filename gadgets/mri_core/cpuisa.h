//
// This is a cross-platform adapataion of the MSDN sample
// https://docs.microsoft.com/en-us/cpp/intrinsics/cpuid-cpuidex
//
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

	const char* CPU_Vendor();
	const char* CPU_Brand();

	bool CPU_supports_SSE3();
	bool CPU_supports_PCLMULQDQ();
	bool CPU_supports_MONITOR();
	bool CPU_supports_SSSE3();
	bool CPU_supports_FMA();
	bool CPU_supports_CMPXCHG16B();
	bool CPU_supports_AVX512POPCNTDQ();
	bool CPU_supports_SSE41();
	bool CPU_supports_SSE42();
	bool CPU_supports_MOVBE();
	bool CPU_supports_POPCNT();
	bool CPU_supports_AES();
	bool CPU_supports_XSAVE();
	bool CPU_supports_OSXSAVE();
	bool CPU_supports_AVX();
	bool CPU_supports_F16C();
	bool CPU_supports_RDRAND();

	bool CPU_supports_MSR();
	bool CPU_supports_CX8();
	bool CPU_supports_SEP();
	bool CPU_supports_CMOV();
	bool CPU_supports_CLFSH();
	bool CPU_supports_MMX();
	bool CPU_supports_FXSR();
	bool CPU_supports_SSE();
	bool CPU_supports_SSE2();

	bool CPU_supports_FSGSBASE();
	bool CPU_supports_BMI1();
	bool CPU_supports_HLE();
	bool CPU_supports_AVX2();
	bool CPU_supports_BMI2();
	bool CPU_supports_ERMS();
	bool CPU_supports_INVPCID();
	bool CPU_supports_RTM();
	bool CPU_supports_AVX512F();
	bool CPU_supports_AVX512DQ();
	bool CPU_supports_RDSEED();
	bool CPU_supports_ADX();
	bool CPU_supports_AVX512IFMA();
	bool CPU_supports_AVX512PF();
	bool CPU_supports_AVX512ER();
	bool CPU_supports_AVX512CD();
	bool CPU_supports_SHA();
	bool CPU_supports_AVX512BW();
	bool CPU_supports_AVX512VL();

	bool CPU_supports_PREFETCHWT1();

	bool CPU_supports_LAHF();
	bool CPU_supports_LZCNT();
	bool CPU_supports_ABM();
	bool CPU_supports_SSE4a();
	bool CPU_supports_XOP();
	bool CPU_supports_TBM();

	bool CPU_supports_SYSCALL();
	bool CPU_supports_MMXEXT();
	bool CPU_supports_RDTSCP();
	bool CPU_supports_3DNOWEXT();
	bool CPU_supports_3DNOW();

#ifdef __cplusplus
}
#endif

