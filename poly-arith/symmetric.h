#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"


#include "fips202.h"

typedef keccak_state xof_state;

#define PQMX_shake128_absorb PQMX_NAMESPACE(PQMX_shake128_absorb)
void PQMX_shake128_absorb(keccak_state *s,
                           const uint8_t seed[PQMX_SYMBYTES],
                           uint8_t x);

#define PQMX_shake256_prf PQMX_NAMESPACE(PQMX_shake256_prf)
void PQMX_shake256_prf(uint8_t *out, size_t outlen, const uint8_t key[PQMX_SYMBYTES], uint8_t nonce);

#define XOF_BLOCKBYTES SHAKE128_RATE

#define hash_h(OUT, IN, INBYTES) sha3_256(OUT, IN, INBYTES)
#define hash_g(OUT, IN, INBYTES) sha3_512(OUT, IN, INBYTES)
#define xof_absorb(STATE, SEED, X) PQMX_shake128_absorb(STATE, SEED, X)
#define xof_squeezeblocks(OUT, OUTBLOCKS, STATE) shake128_squeezeblocks(OUT, OUTBLOCKS, STATE)
#define prf(OUT, OUTBYTES, KEY, NONCE) PQMX_shake256_prf(OUT, OUTBYTES, KEY, NONCE)
#define kdf(OUT, IN, INBYTES) shake256(OUT, PQMX_SSBYTES, IN, INBYTES)


#endif /* SYMMETRIC_H */
