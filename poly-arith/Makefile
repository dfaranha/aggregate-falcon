CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wmissing-prototypes -Wredundant-decls -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer -ggdb -lflint -DDEGREE=64
RM = /bin/rm

SOURCES = poly.c ntt.c reduce.c util.c fips202.c symmetric-shake.c cpucycles.c
HEADERS = params.h poly.h ntt.h reduce.h util.h 

.PHONY: all speed shared clean

all: \
  test_mul

test_mul: $(SOURCES) randombytes.c poly.c ntt.c test_mul.c bench.c
	$(CC) $(CFLAGS) $^ -o test_mul -lm

clean:
	-$(RM) -rf *.gcno *.gcda *.gch *.lcov *.o *.so
	-$(RM) -rf test_mul
