#include <stdint.h>
#include <string.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "symmetric.h"

/* Polynomial defining the cyclotomic ring. */
nmod_poly_t cyclo_poly, crt_poly[2];

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int64_t *r: pointer to output buffer
*              - unsigned int len: requested number of 64-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 64-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int64_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos, j;

  uint64_t t;
  ctr = pos = 0;
  while(ctr < len && pos + 8 <= buflen) {
    t = buf[pos++];
    for(j=8;j<64;j+=8)
      t |= (uint64_t)buf[pos++] << j;
    t &= (1L << 62)-1;

    if(t < PQMX_Q)
      r[ctr++] = t;
  }
  return ctr;
}

/*************************************************
* Name:        poly_uniform
*
* Description: Generate uniform random polynomial 
*
* Arguments:   - const poly *r: pointer to output polynomial
*              - const uint8_t seed[]: pointer to input buffer 
*              (assumed to be uniformly random bytes) of length PQMX_SYMBYTES
*              - uint32_t nonce: 32-bit nonce 
**************************************************/
#define POLY_UNIFORM_NBLOCKS (PQMX_POLYBYTES + XOF_BLOCKBYTES - 1)/(XOF_BLOCKBYTES)
void poly_uniform(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce)
{
  unsigned int i, ctr, off, buflen;
  xof_state state;
  uint8_t rnd[POLY_UNIFORM_NBLOCKS*XOF_BLOCKBYTES+2];
  memset(rnd,0,sizeof(rnd));
  memset(r->coeffs, 0, PQMX_POLYBYTES);

  buflen = POLY_UNIFORM_NBLOCKS*XOF_BLOCKBYTES;

  xof_absorb(&state, seed, nonce);
  xof_squeezeblocks(rnd, POLY_UNIFORM_NBLOCKS, &state);
  
  ctr = rej_uniform(r->coeffs, PQMX_N, rnd, buflen);
  
  while(ctr<PQMX_N) {
    off = buflen % 3;
    for(i = 0; i < off; ++i)
      rnd[i] = rnd[buflen - off + i];
    xof_squeezeblocks(rnd+off, 1, &state);
    buflen = XOF_BLOCKBYTES + off;
    ctr+= rej_uniform(r->coeffs + ctr, PQMX_N - ctr, rnd, buflen);
  }

  poly_reduce(r);
}

/*************************************************
* Name:        poly_ntt
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in bitreversed order
*
* Arguments:   - uint64_t *r: pointer to in/output polynomial
**************************************************/
void poly_ntt4(poly *r)
{
  ntt4(r->coeffs);
  poly_reduce(r);
}

void poly_ntt8(poly *r)
{
  ntt8(r->coeffs);
  poly_reduce(r);
}

void poly_ntt_full(poly *r)
{
  ntt_full(r->coeffs);
  poly_reduce(r);
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT)
*              of a polynomial in place;
*              inputs assumed to be in bitreversed order, output in normal order
*
* Arguments:   - uint64_t *a: pointer to in/output polynomial
**************************************************/
void poly_invntt_tomont(poly *r)
{
  invntt(r->coeffs);
}

/*************************************************
* Name:        poly_basemul_montgomery
*
* Description: Multiplication of two polynomials in NTT domain
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul_montgomery4(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<PQMX_N/8;i++) {
    basemul4(&r->coeffs[8*i], &a->coeffs[8*i], &b->coeffs[8*i], zetas[PQMX_L/2+i]);
    basemul4(&r->coeffs[8*i+4], &a->coeffs[8*i+4], &b->coeffs[8*i+4], -zetas[PQMX_L/2+i]);
  }
}

void poly_basemul_montgomery8(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<PQMX_N/16;i++) {
    basemul8(&r->coeffs[16*i], &a->coeffs[16*i], &b->coeffs[16*i], zetas[PQMX_L/2+i]);
    basemul8(&r->coeffs[16*i+8], &a->coeffs[16*i+8], &b->coeffs[16*i+8], -zetas[PQMX_L/2+i]);
  }
}

/*************************************************
* Name:        poly_pointwise_montgomery
*
* Description: Pointwise multiplication of polynomials in NTT domain
*              representation and multiplication of resulting polynomial
*              by 2^{-32}.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  
  for(i = 0; i < PQMX_N; ++i)
    c->coeffs[i] = montgomery_reduce((int64_t)a->coeffs[i] * b->coeffs[i]);
}

/*************************************************
* Name:        poly_ctrmul
*
* Description: CRT multiplication of two half-degree polynomials.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_crtmul(nmod_poly_t *c, const nmod_poly_t *a, const nmod_poly_t *b) {
	nmod_poly_mulmod(c[0], a[0], b[0], crt_poly[0]);
	nmod_poly_mulmod(c[1], a[1], b[1], crt_poly[1]);
}

/*************************************************
* Name:        poly_tomont
*
* Description: Inplace conversion of all coefficients of a polynomial
*              from normal domain to Montgomery domain
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_tomont(poly *r)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = montgomery_reduce((__int128)r->coeffs[i]*PQMX_MONT2);
}

/*************************************************
* Name:        poly_reduce
*
* Description: Applies Barrett reduction to all coefficients of a polynomial
*              for details of the Barrett reduction see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *r)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = barrett_reduce(r->coeffs[i]);
}

/*************************************************
* Name:        poly_reduce_mont
*
* Description: Applies Montgomery reduction to all coefficients of a polynomial
*              for details of the Montgomery reduction see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_reduce_mont(poly *r)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = montgomery_reduce( (__int128) r->coeffs[i] );
}


/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials; no modular reduction is performed
*
* Arguments: - poly *r: pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials; no modular reduction is performed
*
* Arguments: - poly *r:       pointer to output polynomial+
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}

/*************************************************
* Name:        poly_shift
*
* Description: Inplace Shift polynomial coefficients right by one and negate
*              the new leading coefficient. 
*
* Arguments: - poly *r:       pointer to input polynomial
**************************************************/
void poly_shift(poly *r){
  int64_t tmp = r->coeffs[PQMX_N-1];
  unsigned int i;
  for(i=PQMX_N-1;i>0;i--){
    r->coeffs[i] = r->coeffs[i-1];
  }
  r->coeffs[0] = -1*tmp;
}