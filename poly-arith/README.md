Proof-of-concept code for the polynomial arithmetic experiments accompannying the paper "Aggregating Falcon Signatures With LaBRADOR" submitted to CRYPTO 2024.

Depedencies are the GMP 6.2 and FLINT 2.9 libraries.

For building the code, run `make` inside the source directory. This will build the `test_mul` binary for running tests and benchmarks for various parameters.

By default it run the experiments for subring degree 64, but other parameters can be enabled by replacing `-DDEGREE=64` with `-DDEGREE=128` or `-DDEGREE=256` in the Makefile.
The code implements various polynomial multiplication algorithms for 64-bit modulus. One can then manually scale the latencies to 128-bit modulus.
