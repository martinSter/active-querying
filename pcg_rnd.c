// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for temporal network SIR by Petter Holme (2018)

// this file contains the random number generator, derived from the PCG
// RNG v0.94 http://www.pcg-random.org under the Apache License 2.0
// http://www.apache.org/licenses/LICENSE-2.0

// 32-bit Output, 64-bit State, the PCG-XSH-RR version

#include "run.h"

// declare variable g without defining it (g has been defined elsewhere 
// and here we just declare it so it can be used)
extern GLOBALS g;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t pcg_32 () {
	uint64_t state = g.state;
	uint32_t value, rot;

	g.state = g.state * 6364136223846793005ULL + 1442695040888963407ULL;
	value = ((state >> 18u) ^ state) >> 27u;
	rot = state >> 59u;
	return (value >> rot) | (value << ((- rot) & 31));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t pcg_32_bounded (uint32_t bound) {
	uint32_t threshold = -bound % bound, r;

	for ( ; ; ) {
		r = pcg_32();
		if (r >= threshold) return r % bound;
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint16_t pcg_16 () {
	static unsigned int exist;

	if (exist) {
		exist = 0;
		return g.rmem >> 16;
	}
	exist = 1;
	g.rmem = pcg_32();
	return (uint16_t) g.rmem;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
