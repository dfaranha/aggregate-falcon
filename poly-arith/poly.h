#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

#include <flint/flint.h>
#include <flint/nmod_poly.h>

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct{
  int64_t coeffs[PQMX_N];
} poly;


#define poly_uniform PQMX_NAMESPACE(poly_uniform)
void poly_uniform(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce);

#define poly_ntt4 PQMX_NAMESPACE(poly_ntt4)
void poly_ntt4(poly *r);
#define poly_ntt8 PQMX_NAMESPACE(poly_ntt8)
void poly_ntt8(poly *r);
#define poly_ntt_full PQMX_NAMESPACE(poly_ntt_full)
void poly_ntt_full(poly *r);
#define poly_crtmul PQMX_NAMESPACE(poly_crtmul)
void poly_crtmul(nmod_poly_t *c, const nmod_poly_t *a, const nmod_poly_t *b);
#define poly_invntt_tomont PQMX_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *r);
#define poly_basemul_montgomery4 PQMX_NAMESPACE(poly_basemul_montgomery4)
void poly_basemul_montgomery4(poly *r, const poly *a, const poly *b);
#define poly_basemul_montgomery8 PQMX_NAMESPACE(poly_basemul_montgomery8)
void poly_basemul_montgomery8(poly *r, const poly *a, const poly *b);
#define poly_pointwise_montgomery PQMX_NAMESPACE(poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *r, const poly *a, const poly *b);

#define poly_basemul_acc PQMX_NAMESPACE(poly_basemul_acc)
void poly_basemul_acc(int64_t r[4], const poly *a, const poly *b);

#define poly_tomont PQMX_NAMESPACE(poly_tomont)
void poly_tomont(poly *r);

#define poly_inner_prod_mont PQMX_NAMESPACE(poly_inner_prod_mont)
void poly_inner_prod_mont(poly *r, const poly *a, const poly *b);


#define poly_reduce PQMX_NAMESPACE(poly_reduce)
void poly_reduce(poly *r);
#define poly_reduce_mont PQMX_NAMESPACE(poly_reduce_mont)
void poly_reduce_mont(poly *r);

#define poly_add PQMX_NAMESPACE(poly_add)
void poly_add(poly *r, const poly *a, const poly *b);
#define poly_sub PQMX_NAMESPACE(poly_sub)
void poly_sub(poly *r, const poly *a, const poly *b);

#define poly_shift PQMX_NAMESPACE(poly_shift)
void poly_shift(poly *r);
#endif
