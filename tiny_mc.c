/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 *
 *  To the extent possible under law, the author has dedicated all copyright
 *  and related and neighboring rights to this software to the public domain
 *  worldwide. This software is distributed without any warranty.
 *
 *  See <http://creativecommons.org/publicdomain/zero/1.0/>. */


#define _XOPEN_SOURCE 500  // M_PI
#define N 8

#include <omp.h>
#include "params.h"
#include "wtime.h"

#include <stddef.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <string.h>

/* This is xoshiro512+ 1.0, our generator for floating-point numbers with
 *    increased state size. We suggest to use its upper bits for
 *       floating-point generation, as it is slightly faster than xoshiro512**.
 *          It passes all tests we are aware of except for the lowest three bits,
 *             which might fail linearity tests (and just those), so if low linear
 *                complexity is not considered an issue (as it is usually the case) it
 *                   can be used to generate 64-bit outputs, too.
 *
 *                      We suggest to use a sign test to extract a random Boolean value, and
 *                         right shifts to extract subsets of bits.
 *
 *                            The state must be seeded so that it is not everywhere zero. If you have
 *                               a 64-bit seed, we suggest to seed a splitmix64 generator and use its
 *                                  output to fill s. */

static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}


static uint64_t s[8] = {SEED};

uint64_t next(void) {
	const uint64_t result = s[0] + s[2];

	const uint64_t t = s[1] << 11;

	s[2] ^= s[0];
	s[5] ^= s[1];
	s[1] ^= s[2];
	s[7] ^= s[3];
	s[3] ^= s[4];
	s[4] ^= s[5];
	s[0] ^= s[6];
	s[6] ^= s[7];

	s[6] ^= t;

	s[7] = rotl(s[7], 21);

	return result;
}


/* This is the jump function for the generator. It is equivalent
 *    to 2^256 calls to next(); it can be used to generate 2^256
 *       non-overlapping subsequences for parallel computations. */

void jump(void) {
	static const uint64_t JUMP[] = { 0x33ed89b6e7a353f9, 0x760083d7955323be, 0x2837f2fbb5f22fae, 0x4b8c5674d309511c, 0xb11ac47a7ba28c25, 0xf1be7667092bcc1c, 0x53851efdb6df0aaf, 0x1ebbc8b23eaf25db };

	uint64_t t[sizeof s / sizeof *s];
	memset(t, 0, sizeof t);
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b)
				for(unsigned int w = 0; w < sizeof s / sizeof *s; w++)

					t[w] ^= s[w];
			next();
		}

	memcpy(s, t, sizeof s);	
}


/* This is the long-jump function for the generator. It is equivalent to
 *    2^384 calls to next(); it can be used to generate 2^128 starting points,
 *       from each of which jump() will generate 2^128 non-overlapping
 *          subsequences for parallel distributed computations. */

void long_jump(void) {
	static const uint64_t LONG_JUMP[] = { 0x11467fef8f921d28, 0xa2a819f2e79c8ea8, 0xa8299fc284b3959a, 0xb4d347340ca63ee1, 0x1cb0940bedbff6ce, 0xd956c5c4fa1f8e17, 0x915e38fd4eda93bc, 0x5b3ccdfa5d7daca5 };

	uint64_t t[sizeof s / sizeof *s];
	memset(t, 0, sizeof t);
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b)
				for(unsigned int w = 0; w < sizeof s / sizeof *s; w++)
					t[w] ^= s[w];
			next();
		}

	memcpy(s, t, sizeof s);	
}


/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
" and Nicolas Wolovick";




struct photon {
	float x;
	float y;
	float z;
	float u;
	float v;
	float w;
	float weight;
	float t;
	float xi1;
	float xi2;
	unsigned int shell;
	float heat[SHELLS];
	float heat2[SHELLS];
};

//N photons:
float x[N], y[N], z[N], u[N], v[N], w[N], t[N], xi1[N], xi2[N];
float weight;
unsigned int shell[N];
float heat[SHELLS];
float heat2[SHELLS];

static unsigned int g_seed = SEED;

// global state, heat and heat square in each shell



//################  FUNCIONES PARA FAST_RAND() ###########

// Used to seed the generator. 

void fast_srand(int seed) {
	g_seed = seed;
}

// Compute a pseudorandom integer.
// Output value in range [0, 32767]

int fast_rand(void) {
	g_seed = (214013*g_seed+2531011);
	return (g_seed>>16)&0x7FFF;
}    
//##########################################################


const float albedo = MU_S / (MU_S + MU_A);
const float shells_per_mfp = 1e4 / MICRONS_PER_SHELL / (MU_A + MU_S);

/***
 * Photon
 ***/

void make_photons(void) {
	for (size_t i=0; i<N; ++i) {
		x[i] = .0f;
		y[i] = .0f;
		z[i] = .0f;
		u[i] = .0f;
		v[i] = .0f;
		w[i] = 1.0f;
	}
	weight = 1.0f;
}


void photon_mov(void)
{
	unsigned int vector[N] = {fast_rand(),fast_rand(),fast_rand(),fast_rand(),fast_rand(),fast_rand(),fast_rand(),fast_rand()};
	for(size_t i=0; i<N; i++) {
		t[i] = -logf(vector[i] / (float)32767.0) / (MU_A+MU_S); /* move */
		x[i] += t[i] * u[i];
		y[i] += t[i] * v[i];
		z[i] += t[i] * w[i];

	}
	for(size_t i=0; i<N; i++) {
		shell[i] = sqrtf(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) * shells_per_mfp; /* absorb */
		shell[i] = (shell[i] > SHELLS-1) ? SHELLS-1 : shell[i];
	}

	for(size_t i=0; i<N; i++) {
		heat[shell[i]] += (1.0f - albedo) * weight;
		heat2[shell[i]] += (1.0f - albedo) * (1.0f - albedo) * weight * weight; /* add up squares */
	}
	weight *= albedo;
	for(size_t i=0; i<N; i++) {
		do {
			xi1[i] = 2.0f * fast_rand() / (float)32767.0 - 1.0f;
			xi2[i] = 2.0f * fast_rand() / (float)32767.0 - 1.0f;
			t[i] = xi1[i] * xi1[i] + xi2[i] * xi2[i];
		} while (1.0 < t[i]);
	}
	for(size_t i=0; i<N; i++) {
		u[i] = 2.0f * t[i] - 1.0f;
		v[i] = xi1[i] * sqrtf((1.0f - u[i] * u[i]) / t[i]);
		w[i] = xi2[i] * sqrtf((1.0f - u[i] * u[i]) / t[i]);
	}

}


void photon_move_roullete() {
	for(;;) {
		photon_mov();
		if (weight < 0.001f) { /* roulette */
			if (fast_rand() / (float)32767.0 > 0.1f) //TODO: ver que tanto se puede meter mano aca
				break;
			weight /= 0.1f; //TODO: cambiar division por multiplicacion
		}
	}
}
/***
 * Main matter
 ***/

int main(void)
{
	//photon();
	//long_jump();
	// heading
		printf("# %s\n# %s\n# %s\n", t1, t2, t3);
		printf("# Scattering = %8.3f/cm\n", MU_S);
		printf("# Absorption = %8.3f/cm\n", MU_A);
		printf("# Photons    = %8d\n#\n", PHOTONS);

	// configure RNG
	fast_srand(SEED);
	// start timer
	double start = wtime();
	//simulation
	for (unsigned int i = 0; i < PHOTONS/N; ++i) {
		make_photons();
		for(int k=0; k<72; k++){
			photon_mov();
		}
		photon_move_roullete();
		
	}
	// stop timer
	double end = wtime();
	assert(start <= end);
	double elapsed = end - start;

	printf("# %lf seconds\n", elapsed);
	printf("# %lf K photons per second\n", 1e-3 * PHOTONS / elapsed);

	printf("# Radius\tHeat\n");
	printf("# [microns]\t[W/cm^3]\tError\n");
	float t = 4.0f * M_PI * powf(MICRONS_PER_SHELL, 3.0f) * PHOTONS / 1e12;
	for (unsigned int i = 0; i < SHELLS - 1; ++i) {
	printf("%6.0f\t%12.5f\t%12.5f\n", i * (float)MICRONS_PER_SHELL,
	heat[i] / t / (i * i + i + 1.0 / 3.0),
	sqrt(heat2[i] - heat[i] * heat[i] / PHOTONS) / t / (i * i + i + 1.0f / 3.0f));
	}
	printf("# extra\t%12.5f\n", heat[SHELLS - 1] / PHOTONS);
	
	return 0;
}
