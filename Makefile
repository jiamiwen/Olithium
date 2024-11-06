CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wvla -Wpointer-arith -O3 -march=native -mtune=native
NISTFLAGS += -Wno-unused-result -O3
SOURCES = sign.c packing.c polyvec.c poly.c ntt.c reduce.c rounding.c
SOURCES1 = timetest.c packing.c polyvec.c poly.c ntt.c reduce.c rounding.c
SOURCES2 = signnew.c packing.c polyvec.c poly.c ntt.c reduce.c rounding.c
HEADERS = config.h params.h api.h sign.h packing.h polyvec.h poly.h ntt.h \
  reduce.h rounding.h symmetric.h randombytes.h
KECCAK_SOURCES = $(SOURCES) fips202.c symmetric-shake.c
KECCAK_SOURCES1 = $(SOURCES1) fips202.c symmetric-shake.c
KECCAK_SOURCES2 = $(SOURCES2) fips202.c symmetric-shake.c
KECCAK_HEADERS = $(HEADERS) fips202.h
AES_SOURCES = $(SOURCES) fips202.c aes256ctr.c symmetric-aes.c
AES_HEADERS = $(HEADERS) fips202.h aes256ctr.h


Dilithium_test: timePQC.c rng.c rng.h $(KECCAK_SOURCES1) \
  $(KECCAK_HEADERS)
	$(CC) $(NISTFLAGS) -g -o $@ $< rng.c $(KECCAK_SOURCES1) $(LDFLAGS) -lcrypto

Olithium_test: PQCnew.c rng.c rng.h $(KECCAK_SOURCES2) \
  $(KECCAK_HEADERS)
	$(CC) $(NISTFLAGS) -g -o $@ $< rng.c $(KECCAK_SOURCES2) $(LDFLAGS) -lcrypto

