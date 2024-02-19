#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "randombytes.h"
#include "ntt.h"
#include "poly.h"
#include "assert.h"
#include "bench.h"

#define NTESTS 1
#define N PQMX_N
#define Q PQMX_Q
#define SEEDBYTES PQMX_SYMBYTES

/* Modulus and roots for the two-splitting case. */
#define P 	1092158679064256381
#define P0	267789826476780723
#define P1	824368852587475658

extern nmod_poly_t cyclo_poly, crt_poly[2];
int test(void);
int bench(void);

static void poly_naivemul(poly *c, const poly *a, const poly *b) {
  unsigned int i,j;
  int64_t r[2*N] = {0};

  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      r[i+j] = (r[i+j] + ((__int128)a->coeffs[i]*b->coeffs[j])) % Q;

  for(i = N; i < 2*N; i++)
    r[i-N] = (r[i-N] - r[i]) % Q;

  for(i = 0; i < N; i++)
    c->coeffs[i] = r[i];
}

int test(void) {
  int i, j;
  uint8_t seed[SEEDBYTES];
  uint16_t nonce = 0;
  poly a, b, c, d;

  randombytes(seed, sizeof(seed));
  for(i = 0; i < NTESTS; ++i) {
    poly_uniform(&a, seed, nonce++);
    poly_uniform(&b, seed, nonce++);

    c = a;
#if DEGREE < 256
    poly_ntt4(&c);
    for(j = 0; j < N; ++j)
      c.coeffs[j] = ((__int128)c.coeffs[j]*PQMX_MONT3) % Q;
    poly_invntt_tomont(&c);
    for(j = 0; j < N; ++j) {
      if((c.coeffs[j] - a.coeffs[j]) % Q)
        fprintf(stderr, "ERROR in ntt/invntt: c[%d] = %ld != %ld\n", j, c.coeffs[j]%Q, a.coeffs[j]);
    }

    poly_naivemul(&c, &a, &b);
    poly_ntt4(&a);
    poly_ntt4(&b);
    poly_basemul_montgomery4(&d, &a, &b);
    poly_invntt_tomont(&d);
#else
    poly_ntt8(&c);
    for(j = 0; j < N; ++j)
      c.coeffs[j] = ((__int128)c.coeffs[j]*PQMX_MONT3) % Q;
    poly_invntt_tomont(&c);
    for(j = 0; j < N; ++j) {
      if((c.coeffs[j] - a.coeffs[j]) % Q)
        fprintf(stderr, "ERROR in ntt/invntt: c[%d] = %ld != %ld\n", j, c.coeffs[j]%Q, a.coeffs[j]);
    }

    poly_naivemul(&c, &a, &b);
    poly_ntt8(&a);
    poly_ntt8(&b);
    poly_basemul_montgomery8(&d, &a, &b);
    poly_invntt_tomont(&d);
#endif

    for(j = 0; j < N; ++j) {
      if((d.coeffs[j] - c.coeffs[j]) % Q)
        fprintf(stderr, "ERROR in multiplication: d[%d] = %ld != %ld\n", j, d.coeffs[j], c.coeffs[j]);
    }
  }

  return 0;
}

int bench(void) {
  uint8_t seed[SEEDBYTES];
  uint16_t nonce = 0;
  poly a, b, c, d;
  nmod_poly_t crt_a[2], crt_b[2], crt_c[2];
  uint64_t coeff;

  for (size_t i = 0; i < 2; i++) {
    nmod_poly_init(crt_a[i], P);
    nmod_poly_init(crt_b[i], P);
    nmod_poly_init(crt_c[i], P);
  }

  randombytes(seed, sizeof(seed));
  poly_uniform(&a, seed, nonce++);
  poly_uniform(&b, seed, nonce++);
  for (int i = 0; i < N/2; i++) {
	for (int j = 0; j < 2; j++) {
      randombytes((uint8_t *)&coeff, sizeof(coeff));
      nmod_poly_set_coeff_ui(crt_a[j], i, coeff % P);
      randombytes((uint8_t *)&coeff, sizeof(coeff));
      nmod_poly_set_coeff_ui(crt_b[j], i, coeff % P);
      randombytes((uint8_t *)&coeff, sizeof(coeff));
      nmod_poly_set_coeff_ui(crt_c[j], i, coeff % P);
	}
  }

  BENCH_SMALL("field mul:               ", field_mul(&c.coeffs[0], &a.coeffs[0], &b.coeffs[0]));
  BENCH_SMALL("almost-full up-to-4 ntt: ", poly_ntt4(&a));
  BENCH_SMALL("almost-full up-to-8 ntt: ", poly_ntt8(&a));
  BENCH_SMALL("fully-split ntt:         ", poly_ntt_full(&a));
  BENCH_SMALL("inv ntt:                 ", poly_invntt_tomont(&d);)
  BENCH_SMALL("naive schoolbook mul:    ", poly_naivemul(&c, &a, &b));
  BENCH_SMALL("two-split CRT mul:       ", poly_crtmul(crt_c, crt_a, crt_b));
  BENCH_SMALL("almost-full up-to-4 mul: ", poly_basemul_montgomery4(&d, &a, &b););
  BENCH_SMALL("almost-full up-to-8 mul: ", poly_basemul_montgomery8(&d, &a, &b););
  BENCH_SMALL("fully-split mul:         ", poly_pointwise_montgomery(&c, &a, &b));

  for (size_t i = 0; i < 2; i++) {
    nmod_poly_clear(crt_a[i]);
    nmod_poly_clear(crt_b[i]);
    nmod_poly_clear(crt_c[i]);
  }

  return 0;
}

int main(void) {
	nmod_poly_init(cyclo_poly, P);
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(crt_poly[i], P);
	}

	// Initialize polynomial as x^N + 1. */
	nmod_poly_set_coeff_ui(cyclo_poly, N, 1);
	nmod_poly_set_coeff_ui(cyclo_poly, 0, 1);

	// Initialize two factors of the polynomial for CRT representation.
	nmod_poly_set_coeff_ui(crt_poly[0], N/2, 1);
	nmod_poly_set_coeff_ui(crt_poly[0], 0, P0);
	nmod_poly_set_coeff_ui(crt_poly[1], N/2, 1);
	nmod_poly_set_coeff_ui(crt_poly[1], 0, P1);

	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(crt_poly[i]);
	}
	nmod_poly_clear(cyclo_poly);

	test();
	bench();
}