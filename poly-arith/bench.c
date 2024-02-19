#include <stdio.h>
#include "bench.h"
#include "cpucycles.h"

// Private Definitions

//Stores the time measured before the execution of the benchmark.
static long long before;

//Stores the time measured after the execution of the benchmark.
static long long after;

//Stores the sum of timings for the current benchmark.
static long long total;

// Public Definitions

void bench_reset() {
	total = 0;
}

void bench_before() {
	before = cpucycles();
}

void bench_after() {
	long long result;
	after = cpucycles();
	result = (after - before);
	total += result;
}

void bench_compute(int benches) {
	total = total / benches;
}

void bench_print() {
	printf("%lld cycles\n", total);
	printf("\n");
}

unsigned long long bench_total() {
	return total;
}
