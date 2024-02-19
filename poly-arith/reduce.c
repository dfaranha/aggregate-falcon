#include <stdint.h>
#include "params.h"
#include "reduce.h"

/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 128-bit integer a, computes
*              64-bit integer congruent to a * R^-1 mod q, where R=2^64
*
* Arguments:   - __int128 a: input integer to be reduced;
*                           has to be in {-q2^63,...,q2^63-1}
*
* Returns:     integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
**************************************************/
int64_t montgomery_reduce(__int128 a)
{
  int64_t t;

  t = (int64_t)a*PQMX_QINV;
  t = (a - (__int128)t*PQMX_Q) >> 64;
  return t;
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 64-bit integer a, computes
*              centered representative congruent to a mod q in {-(q-1)/2,...,(q-1)/2}
*
* Arguments:   - int64_t a: input integer to be reduced
*
* Returns:     integer in {-(q-1)/2,...,(q-1)/2} congruent to a modulo q.
**************************************************/
int64_t barrett_reduce(int64_t a) {
  int64_t t;
  const int64_t v = ( ((__int128)1 <<124 ) + PQMX_Q/2)/PQMX_Q;
  t  = ((__int128)v*a + ((__int128)1 << 123) ) >> 124;
  t *= PQMX_Q;
  return a - t;
}