#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define zetas PQMX_NAMESPACE(zetas)
extern const int64_t zetas[PQMX_N];

#define ntt4 PQMX_NAMESPACE(ntt4)
void ntt4(int64_t poly[PQMX_N]);

#define ntt8 PQMX_NAMESPACE(ntt8)
void ntt8(int64_t poly[PQMX_N]);

#define ntt_full PQMX_NAMESPACE(ntt_full)
void ntt_full(int64_t poly[PQMX_N]);

#define invntt PQMX_NAMESPACE(invntt)
void invntt(int64_t poly[PQMX_N]);

#define basemul4 PQMX_NAMESPACE(basemul4)
void basemul4(int64_t r[4], const int64_t a[4], const int64_t b[4], int64_t zeta);

#define basemul8 PQMX_NAMESPACE(basemul8)
void basemul8(int64_t r[8], const int64_t a[8], const int64_t b[8], int64_t zeta);

#define scalar_field_mul PQMX_NAMESPACE(scalar_field_mul)
void scalar_field_mul(int64_t r[4], const int64_t a, const int64_t b[4]);

#define field_mul PQMX_NAMESPACE(field_mul)
void field_mul(int64_t *c, int64_t *a, int64_t *b);

#define fqmul PQMX_NAMESPACE(fqmul)
int64_t fqmul(int64_t a, int64_t b);

#endif
