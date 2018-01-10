#ifndef __asiosys__
	#define __asiosys__

	#ifdef _WIN32
		#undef MAC
		#undef _MAC
		#define PPC 0
		#define WINDOWS 1
		#define SGI 0
		#define SUN 0
		#define LINUX 0
		#define BEOS 0

		#define NATIVE_INT64 0
		#define IEEE754_64FLOAT 1

#endif

#endif
