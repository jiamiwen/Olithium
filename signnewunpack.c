#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"

#include<stdio.h>
#include<time.h>
/*************************************************
* Name:        crypto_sign_keypair
*
* Description: Generates public and private key.
*
* Arguments:   - uint8_t *pk: pointer to output public key (allocated
*                             array of CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key (allocated
*                             array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk) {
  uint8_t seedbuf[3*SEEDBYTES];
  uint8_t tr[CRHBYTES];
  const uint8_t *rho, *rhoprime, *key;
  polyvecl mat[K];
  polyvecl s1, s1hat;
  polyveck s2, t1, t0;

  /* Get randomness for rho, rhoprime and key */
  randombytes(seedbuf, SEEDBYTES);
  shake256(seedbuf, 3*SEEDBYTES, seedbuf, SEEDBYTES);
  rho = seedbuf;
  rhoprime = seedbuf + SEEDBYTES;
  key = seedbuf + 2*SEEDBYTES;

  /* Expand matrix */
  polyvec_matrix_expand(mat, rho);

  /* Sample short vectors s1 and s2 */
  polyvecl_uniform_eta(&s1, rhoprime, 0);
  polyveck_uniform_eta(&s2, rhoprime, L);

  /* Matrix-vector multiplication */
  s1hat = s1;
  polyvecl_ntt(&s1hat);
  polyvec_matrix_pointwise_montgomery(&t1, mat, &s1hat);
  polyveck_reduce(&t1);
  polyveck_invntt_tomont(&t1);

  /* Add error vector s2 */
  polyveck_add(&t1, &t1, &s2);

  /* Extract t1 and write public key */
  polyveck_caddq(&t1);
  polyveck_power2round(&t1, &t0, &t1);
  pack_pk(pk, rho, &t1);

  /* Compute CRH(rho, t1) and write secret key */
  crh(tr, pk, CRYPTO_PUBLICKEYBYTES);
  pack_sk(sk, rho, tr, key, &t0, &s1, &s2);

  return 0;
}

void xorArrays(unsigned char* arr1, unsigned char* arr2, int len) {
	for (int i = 0; i < len; i++) {
		arr1[i] ^= arr2[i];
	}
}
void xorArrays1(unsigned char* arr0, unsigned char* arr1, unsigned char* arr2, int len) {
	for (int i = 0; i < len; i++) {
		arr0[i] = arr1[i] ^ arr2[i];
	}
}

void printBstr(char *S, unsigned char *A, unsigned long long l)
{
	unsigned long long  i;

	printf("%s", S);

	for ( i=0; i<l; i++ )
		printf("%02X", A[i]);

	if ( L == 0 )
		printf("00");

	printf("\n");
}
/*************************************************
* Name:        crypto_sign_signature
*
* Description: Computes signature.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk)
{
  unsigned int n;
  uint8_t seedbuf[2*SEEDBYTES + 3*CRHBYTES];
  uint8_t c[SEEDBYTES];
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;//表示小k
  polyvecl mat[K], s1, y, z;//l纬的向量
  polyveck t0, s2, w1, w0, h;//k纬的向量  这里h是k纬向量
  poly cp;
  keccak_state state;

  const int ee=30;
  polyvecl ytemp[ee];
  uint8_t ctemp[SEEDBYTES*ee];
  polyveck w0temp[ee];
  polyveck w1temp[ee];

  rho = seedbuf;//占32个字节
  tr = rho + SEEDBYTES;// 2019年版本 占48个字节
  key = tr + CRHBYTES;//32个字节
  mu = key + SEEDBYTES;//占48个字节
  rhoprime = mu + CRHBYTES;//占48个字节
  unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);//对私钥解包

  polyvec_matrix_expand(mat, rho);//ExpandA(rho)

  /* Compute CRH(tr, msg) */
  shake256_init(&state);
  shake256_absorb(&state, tr, CRHBYTES);
  //shake256_absorb(&state, m, mlen); 修改处
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else
  crh(rhoprime, key, SEEDBYTES + CRHBYTES);//传入一个地址 因为key和mu的存储位置相邻 CRH(key||mu)
#endif
  polyvecl_ntt(&s1); //Pre-compute
  polyveck_ntt(&s2);//Pre-compute
  polyveck_ntt(&t0);//Pre-compute
  
//rej:
 /* Sample intermediate vector y */
 for(int t=0;t<ee;t++){
  polyvecl_uniform_gamma1(&ytemp[t], rhoprime, nonce++);//y=ExpandMask(rhoprime,k)
  z = ytemp[t];
  polyvecl_ntt(&z);
  /* Matrix-vector multiplication */
  polyvec_matrix_pointwise_montgomery(&w1temp[t], mat, &z);//w1=Ay
  polyveck_reduce(&w1temp[t]);//[-6283009,6283007].
  polyveck_invntt_tomont(&w1temp[t]);//invntt 如果没有 pointwise_montgomery 最后的结果就会在蒙哥马利域中 

  /* Decompose w and call the random oracle */
  polyveck_caddq(&w1temp[t]); //For all coefficients of polynomials in vector of length K add q if coefficient is negative.
  polyveck_decompose(&w1temp[t], &w0temp[t], &w1temp[t]); //HighBits(w1) LowBits(w1) input 这里w0=r0
  polyveck_pack_w1(sig, &w1temp[t]);//暂时将w1打包进签名 哈希的时候需要用打包格式

  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, sig, K*POLYW1_PACKEDBYTES);//sig为w1
  shake256_finalize(&state);
  shake256_squeeze(ctemp+(SEEDBYTES*t) ,SEEDBYTES, &state);//sig=CRH(mu||w1) 也就是c波浪
  }
  //poly_challenge(&cp, c);//c=SampleInBall(c波浪)
  /* Expand matrix and transform vectors */
  //poly_ntt(&cp);
  int ww;
  clock_t start_t, finish_t;
  double total_t;
  int i;
  start_t = clock();
  for(ww=0;ww<2000;ww++){
  i=-1;
  shake256_init(&state);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(sig, SEEDBYTES, &state);//sig=H(M) 也就是k波浪
rej:
  i++;
  xorArrays1(c,ctemp+(SEEDBYTES*i) ,sig,SEEDBYTES);
  poly_challenge(&cp, c);//k=SampleInBall(k波浪)
  //printBstr("c:", ctemp+(SEEDBYTES*i), SEEDBYTES);
  poly_ntt(&cp);
  /* Compute z, reject if it reveals secret */
  polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
  polyvecl_invntt_tomont(&z);
  polyvecl_add(&z, &z, &ytemp[i]);
  polyvecl_reduce(&z);
  if(polyvecl_chknorm(&z, GAMMA1 - BETA))//z的无穷范数大于GAMMA1 - BETA1
    goto rej;

  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
  polyveck_invntt_tomont(&h);
  polyveck_sub(&w0, &w0temp[i], &h);
  polyveck_reduce(&w0);
  if(polyveck_chknorm(&w0, GAMMA2 - BETA))
    goto rej;

  /* Compute hints for w1 */
  polyveck_pointwise_poly_montgomery(&h, &cp, &t0);
  polyveck_invntt_tomont(&h);
  polyveck_reduce(&h);
  if(polyveck_chknorm(&h, GAMMA2))
    goto rej;

  polyveck_add(&w0, &w0, &h);
  polyveck_caddq(&w0); 
  n = polyveck_make_hint(&h, &w0, &w1temp[i]);
  if(n > OMEGA)
    goto rej;
  }
  /* Write signature */
  pack_sig(sig,c , &z, &h);//签名写入c、z、h  sig暂存了c
  *siglen = CRYPTO_BYTES;
  finish_t = clock();
  total_t = (double)(finish_t - start_t) / CLOCKS_PER_SEC;//将时间转换为秒
	printf("CPU 占用的总时间：%f\n", total_t/ww);
  printf("拒绝次数为：%d\n", i+1);
  return 0;
}

/*************************************************
* Name:        crypto_sign
*
* Description: Compute signed message.
*
* Arguments:   - uint8_t *sm: pointer to output signed message (allocated
*                             array with CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smlen: pointer to output length of signed
*                               message
*              - const uint8_t *m: pointer to message to be signed
*              - size_t mlen: length of message
*              - const uint8_t *sk: pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk)
{
  size_t i;

  for(i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, sk);
  *smlen += mlen;
  return 0;
}

/*************************************************
* Name:        crypto_sign_verify
*
* Description: Verifies signature.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_verify(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk)
{
  unsigned int i;
  uint8_t buf[K*POLYW1_PACKEDBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[SEEDBYTES];
  uint8_t c2[SEEDBYTES];
  poly cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h;
  keccak_state state;

  uint8_t k[SEEDBYTES];

  if(siglen != CRYPTO_BYTES)
    return -1;

  unpack_pk(rho, &t1, pk);
  if(unpack_sig(c, &z, &h, sig))
    return -1;
  if(polyvecl_chknorm(&z, GAMMA1 - BETA))
    return -1;

  /* Compute CRH(CRH(rho, t1), msg) */
  crh(mu, pk, CRYPTO_PUBLICKEYBYTES);
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  //shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  poly_challenge(&cp, c);
  polyvec_matrix_expand(mat, rho);

  polyvecl_ntt(&z);
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z);

  poly_ntt(&cp);
  polyveck_shiftl(&t1);
  polyveck_ntt(&t1);
  polyveck_pointwise_poly_montgomery(&t1, &cp, &t1);

  polyveck_sub(&w1, &w1, &t1);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont(&w1);

  /* Reconstruct w1 */
  polyveck_caddq(&w1);
  polyveck_use_hint(&w1, &w1, &h);
  polyveck_pack_w1(buf, &w1);

  /* Call random oracle and verify challenge */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, buf, K*POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(c2, SEEDBYTES, &state);
  
  shake256_init(&state);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(k, SEEDBYTES, &state);//sig=H(M) 也就是k波浪

  xorArrays(c2,k,SEEDBYTES);

  for(i = 0; i < SEEDBYTES; ++i)
    if(c[i] != c2[i])
      return -1;

  return 0;
}

/*************************************************
* Name:        crypto_sign_open
*
* Description: Verify signed message.
*
* Arguments:   - uint8_t *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mlen: pointer to output length of message
*              - const uint8_t *sm: pointer to signed message
*              - size_t smlen: length of signed message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_open(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *pk)
{
  size_t i;

  if(smlen < CRYPTO_BYTES)
    goto badsig;

  *mlen = smlen - CRYPTO_BYTES;
  if(crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk))
    goto badsig;
  else {
    /* All good, copy msg, return 0 */
    for(i = 0; i < *mlen; ++i)
      m[i] = sm[CRYPTO_BYTES + i];
    return 0;
  }

badsig:
  /* Signature verification failed */
  *mlen = -1;
  for(i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}
