#include "util.h"
#include <math.h>
#include <string.h>
#include "randombytes.h"
#include "reduce.h"

/*************************************************
* Name:        poly_norm
*
* Description: Calculate l_p norm of a polynomial
*
* Arguments:   - const poly *a: pointer to input polynomial
*              - int l : norm degree. 0 for infinity norm.  
**************************************************/
double poly_norm(const poly *a, int l)
{
    double norm = 0.0;
    int64_t i;
    double j;
    for(i=0;i<PQMX_N;i++){
        j = fabs((double) a->coeffs[i]);
        if(l == 0)
            norm = max(norm, j);
        else
            norm += pow(j, l*1.0);
    }
    if(l==0) return norm;
    return pow(norm, 1.0/l);
}

/*************************************************
* Name:        polyvec_norm
*
* Description: Calculate l_p norm of a polynomial vector
*
* Arguments:   - const poly *a: pointer to input polynomial
*              - int l : norm degree. 0 for infinity norm. 
*              - int vlen: vector dimension 
**************************************************/
double polyvec_norm(const poly *a, int l, int len)
{
    double norm = 0.0, tmp;
    int64_t i;

    for(i=0;i<len;i++){
        tmp = poly_norm(&a[i], l);
        if(l == 0)
            norm = max(norm, tmp);
        else
            norm += pow(tmp, l*1.0);
    }
    if(l==0) return norm;
    return pow(norm, 1.0/l);
}

/*************************************************
* Name:        secure_random
*
* Description: Generate secure random number between 1 and bound
*
* Arguments:   - int j: bound 
**************************************************/
int secure_random(int j)
{
	uint8_t a[4];
	randombytes(a, 4);
	uint32_t * _a = (uint32_t *)a;
	return *_a % j;
}

/*************************************************
* Name:        shuffle_array
*
* Description: Shuffle array using naive shuffle approach
*
* Arguments:   - int64_t *arr: pointer to input array
*              - int len: input array length
**************************************************/
void shuffle_array(int64_t *arr, int len)
{
    int i,j;
    int64_t tmp;
    for(i=len-1; i > 0; i--){
        j = secure_random(i);
        tmp = arr[j];
        arr[j] = arr[i];
        arr[i] = tmp;
    }
}

/*************************************************
* Name:        padding
*
* Description: A simple padding mechanism. padded-data =  \x01 \x01 *\xFF \x01 data \x00
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - const uint8_t *buf: pointer to input buffer
*              - size_t buflen: buffer length
**************************************************/
void padding(uint8_t out[PQMX_INDCPA_MSGBYTES], const uint8_t *buf, size_t buflen)
{
    memset(out, 0xFF, PQMX_INDCPA_MSGBYTES);
    out[0] = 0x01;
    out[1] = 0x01;
    out[PQMX_INDCPA_MSGBYTES-1] = 0x00;
    out[PQMX_INDCPA_MSGBYTES - buflen - 1] = 0x01;
    memcpy(out+PQMX_INDCPA_MSGBYTES-buflen, buf, buflen);
}

/*************************************************
* Name:        remove_padding
*
* Description: Unpadding
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - const uint8_t *in: pointer to input buffer
**************************************************/
void remove_padding(uint8_t *out, const uint8_t in[PQMX_INDCPA_MSGBYTES])
{
    memset(out, 0, 32); 
    int i=2; 
    while(in[i] == 0xFF){
        i++;
    }
    memcpy(out, in+i+1, PQMX_INDCPA_MSGBYTES - i-1);    
}