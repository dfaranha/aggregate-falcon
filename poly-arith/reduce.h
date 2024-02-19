#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "params.h"

#define montgomery_reduce PQMX_NAMESPACE(montgomery_reduce)
int64_t montgomery_reduce(__int128 a);

#define barrett_reduce PQMX_NAMESPACE(barrett_reduce)
int64_t barrett_reduce(int64_t a);

#endif
