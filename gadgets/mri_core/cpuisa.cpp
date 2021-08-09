#include <string.h>

#ifdef _MSC_VER
// SIMD intrinsics for Windows
#include <intrin.h>
#define cpuid(info, level) __cpuid(info, level)
#define cpuidex(info, leaf, subleaf) __cpuidex(info, leaf, subleaf)
#else
// SIMD intrinsics for GCC
#include <x86intrin.h>
#include <cpuid.h>
#define cpuid(info, level) __cpuid(level, info[0], info[1], info[2], info[3]);
#define cpuidex(info, leaf, subleaf) __cpuid_count(leaf, subleaf, info[0], info[1], info[2], info[3]);
#endif

#include "cpuisa.h"

namespace {
	class BitSet
	{
		int bits;

	public:
		BitSet()
			: bits(0) {}

		void operator = (const int v)
		{
			bits = v;
		}

		bool operator[] (int i) const
		{
			return (bits >> i) & 0x00000001;
		}
	};

	class DataSet
	{
		int* data;

	public:
		DataSet(const int highestId)
		{
			data = (int*)calloc((highestId + 1) << 2, sizeof(int));
		}

		~DataSet()
		{
			free(data);
		}

		int* operator[] (unsigned int i) const
		{
			return (data + (i << 2));
		}
	};

	class InstructionSet
	{
	public:
		InstructionSet();

		char vendor_[0x20];
		char brand_[0x40];
		bool isIntel_;
		bool isAMD_;
		BitSet f_1_ECX_;
		BitSet f_1_EDX_;
		BitSet f_7_EBX_;
		BitSet f_7_ECX_;
		BitSet f_81_ECX_;
		BitSet f_81_EDX_;
	};

	InstructionSet::InstructionSet()
		: isIntel_{ false },
		isAMD_{ false },
		f_1_ECX_{},
		f_1_EDX_{},
		f_7_EBX_{},
		f_7_ECX_{},
		f_81_ECX_{},
		f_81_EDX_{}
	{
		int cpui[4] = { -1 };

		// Calling __cpuid with 0x0 as the function_id argument  
		// gets the number of the highest valid function ID.  
		cpuid(cpui, 0);
		unsigned int nIds_ = cpui[0];

		DataSet data_(nIds_);
		for (unsigned int i = 0; i <= nIds_; ++i)
		{
			cpuidex(data_[i], i, 0);
		}

		// Capture vendor string  
		memset(vendor_, 0, sizeof(vendor_));
		*reinterpret_cast<int*>(vendor_) = data_[0][1];
		*reinterpret_cast<int*>(vendor_ + 4) = data_[0][3];
		*reinterpret_cast<int*>(vendor_ + 8) = data_[0][2];

		if (strncmp(vendor_, "GenuineIntel", 0x20) == 0)
		{
			isIntel_ = true;
		}
		else if (strncmp(vendor_, "AuthenticAMD", 0x20) == 0)
		{
			isAMD_ = true;
		}

		// load bitset with flags for function 0x00000001  
		if (nIds_ >= 1)
		{
			f_1_ECX_ = data_[1][2];
			f_1_EDX_ = data_[1][3];
		}

		// load bitset with flags for function 0x00000007  
		if (nIds_ >= 7)
		{
			f_7_EBX_ = data_[7][1];
			f_7_ECX_ = data_[7][2];
		}

		// Calling __cpuid with 0x80000000 as the function_id argument  
		// gets the number of the highest valid extended ID.  
		cpuid(cpui, 0x80000000);
		unsigned int nExIds_ = cpui[0];

		DataSet extdata_(nExIds_ - 0x80000000);
		for (unsigned int i = 0x80000000; i <= nExIds_; ++i)
		{
			cpuidex(extdata_[i - 0x80000000], i, 0);
		}

		// load bitset with flags for function 0x80000001  
		if (nExIds_ >= 0x80000001)
		{
			f_81_ECX_ = extdata_[1][2];
			f_81_EDX_ = extdata_[1][3];
		}

		memset(brand_, 0, sizeof(brand_));

		// Interpret CPU brand string if reported  
		if (nExIds_ >= 0x80000004)
		{
			memcpy(brand_, extdata_[2], sizeof(int) * 12);
		}
	};

	static const InstructionSet CPU_Rep;
}

const char* CPU_Vendor() { return reinterpret_cast<const char*>(CPU_Rep.vendor_); }
const char* CPU_Brand() { return reinterpret_cast<const char*>(CPU_Rep.brand_); }

bool CPU_supports_SSE3() { return CPU_Rep.f_1_ECX_[0]; }
bool CPU_supports_PCLMULQDQ() { return CPU_Rep.f_1_ECX_[1]; }
bool CPU_supports_MONITOR() { return CPU_Rep.f_1_ECX_[3]; }
bool CPU_supports_SSSE3() { return CPU_Rep.f_1_ECX_[9]; }
bool CPU_supports_FMA() { return CPU_Rep.f_1_ECX_[12]; }
bool CPU_supports_CMPXCHG16B() { return CPU_Rep.f_1_ECX_[13]; }
bool CPU_supports_AVX512POPCNTDQ() { return CPU_Rep.f_1_ECX_[14]; }
bool CPU_supports_SSE41() { return CPU_Rep.f_1_ECX_[19]; }
bool CPU_supports_SSE42() { return CPU_Rep.f_1_ECX_[20]; }
bool CPU_supports_MOVBE() { return CPU_Rep.f_1_ECX_[22]; }
bool CPU_supports_POPCNT() { return CPU_Rep.f_1_ECX_[23]; }
bool CPU_supports_AES() { return CPU_Rep.f_1_ECX_[25]; }
bool CPU_supports_XSAVE() { return CPU_Rep.f_1_ECX_[26]; }
bool CPU_supports_OSXSAVE() { return CPU_Rep.f_1_ECX_[27]; }
bool CPU_supports_AVX() { return CPU_Rep.f_1_ECX_[28]; }
bool CPU_supports_F16C() { return CPU_Rep.f_1_ECX_[29]; }
bool CPU_supports_RDRAND() { return CPU_Rep.f_1_ECX_[30]; }

bool CPU_supports_MSR() { return CPU_Rep.f_1_EDX_[5]; }
bool CPU_supports_CX8() { return CPU_Rep.f_1_EDX_[8]; }
bool CPU_supports_SEP() { return CPU_Rep.f_1_EDX_[11]; }
bool CPU_supports_CMOV() { return CPU_Rep.f_1_EDX_[15]; }
bool CPU_supports_CLFSH() { return CPU_Rep.f_1_EDX_[19]; }
bool CPU_supports_MMX() { return CPU_Rep.f_1_EDX_[23]; }
bool CPU_supports_FXSR() { return CPU_Rep.f_1_EDX_[24]; }
bool CPU_supports_SSE() { return CPU_Rep.f_1_EDX_[25]; }
bool CPU_supports_SSE2() { return CPU_Rep.f_1_EDX_[26]; }

bool CPU_supports_FSGSBASE() { return CPU_Rep.f_7_EBX_[0]; }
bool CPU_supports_BMI1() { return CPU_Rep.f_7_EBX_[3]; }
bool CPU_supports_HLE() { return CPU_Rep.isIntel_ && CPU_Rep.f_7_EBX_[4]; }
bool CPU_supports_AVX2() { return CPU_Rep.f_7_EBX_[5]; }
bool CPU_supports_BMI2() { return CPU_Rep.f_7_EBX_[8]; }
bool CPU_supports_ERMS() { return CPU_Rep.f_7_EBX_[9]; }
bool CPU_supports_INVPCID() { return CPU_Rep.f_7_EBX_[10]; }
bool CPU_supports_RTM() { return CPU_Rep.isIntel_ && CPU_Rep.f_7_EBX_[11]; }
bool CPU_supports_AVX512F() { return CPU_Rep.f_7_EBX_[16]; }
bool CPU_supports_AVX512DQ() { return CPU_Rep.f_7_EBX_[17]; }
bool CPU_supports_RDSEED() { return CPU_Rep.f_7_EBX_[18]; }
bool CPU_supports_ADX() { return CPU_Rep.f_7_EBX_[19]; }
bool CPU_supports_AVX512IFMA() { return CPU_Rep.f_7_EBX_[17]; }
bool CPU_supports_AVX512PF() { return CPU_Rep.f_7_EBX_[26]; }
bool CPU_supports_AVX512ER() { return CPU_Rep.f_7_EBX_[27]; }
bool CPU_supports_AVX512CD() { return CPU_Rep.f_7_EBX_[28]; }
bool CPU_supports_SHA() { return CPU_Rep.f_7_EBX_[29]; }
bool CPU_supports_AVX512BW() { return CPU_Rep.f_7_EBX_[30]; }
bool CPU_supports_AVX512VL() { return CPU_Rep.f_7_EBX_[31]; }

bool CPU_supports_PREFETCHWT1() { return CPU_Rep.f_7_ECX_[0]; }

bool CPU_supports_LAHF() { return CPU_Rep.f_81_ECX_[0]; }
bool CPU_supports_LZCNT() { return CPU_Rep.isIntel_ && CPU_Rep.f_81_ECX_[5]; }
bool CPU_supports_ABM() { return CPU_Rep.isAMD_ && CPU_Rep.f_81_ECX_[5]; }
bool CPU_supports_SSE4a() { return CPU_Rep.isAMD_ && CPU_Rep.f_81_ECX_[6]; }
bool CPU_supports_XOP() { return CPU_Rep.isAMD_ && CPU_Rep.f_81_ECX_[11]; }
bool CPU_supports_TBM() { return CPU_Rep.isAMD_ && CPU_Rep.f_81_ECX_[21]; }

bool CPU_supports_SYSCALL() { return CPU_Rep.isIntel_ && CPU_Rep.f_81_EDX_[11]; }
bool CPU_supports_MMXEXT() { return CPU_Rep.isAMD_ && CPU_Rep.f_81_EDX_[22]; }
bool CPU_supports_RDTSCP() { return CPU_Rep.isIntel_ && CPU_Rep.f_81_EDX_[27]; }
bool CPU_supports_3DNOWEXT() { return CPU_Rep.isAMD_ && CPU_Rep.f_81_EDX_[30]; }
bool CPU_supports_3DNOW() { return CPU_Rep.isAMD_ && CPU_Rep.f_81_EDX_[31]; }


