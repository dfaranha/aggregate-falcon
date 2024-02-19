#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "params.h"
#include "symmetric.h"
#include "fips202.h"

/*************************************************
* Name:        PQMX_shake128_absorb
*
* Description: Absorb step of the SHAKE128 specialized for the PQMX context.
*
* Arguments:   - keccak_state *state: pointer to (uninitialized) output Keccak state
*              - const uint8_t *seed: pointer to PQMX_SYMBYTES input to be absorbed into state
*              - uint8_t x: additional byte of input
**************************************************/
void PQMX_shake128_absorb(keccak_state *state,
                           const uint8_t seed[PQMX_SYMBYTES],
                           uint8_t x)
{
  uint8_t extseed[PQMX_SYMBYTES+1];

  memcpy(extseed, seed, PQMX_SYMBYTES);
  extseed[PQMX_SYMBYTES+0] = x;
  //extseed[PQMX_SYMBYTES+1] = y;

  shake128_absorb_once(state, extseed, sizeof(extseed));
}

/*************************************************
* Name:        PQMX_shake256_prf
*
* Description: Usage of SHAKE256 as a PRF, concatenates secret and public input
*              and then generates outlen bytes of SHAKE256 output
*
* Arguments:   - uint8_t *out: pointer to output
*              - size_t outlen: number of requested output bytes
*              - const uint8_t *key: pointer to the key (of length PQMX_SYMBYTES)
*              - uint8_t nonce: single-byte nonce (public PRF input)
**************************************************/
void PQMX_shake256_prf(uint8_t *out, size_t outlen, const uint8_t key[PQMX_SYMBYTES], uint8_t nonce)
{
  uint8_t extkey[PQMX_SYMBYTES+1];

  memcpy(extkey, key, PQMX_SYMBYTES);
  extkey[PQMX_SYMBYTES] = nonce;

  shake256(out, outlen, extkey, sizeof(extkey));
}
