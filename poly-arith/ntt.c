#include <stdint.h>
#include "params.h"
#include "ntt.h"
#include "reduce.h"

// generate with the following SageMath code:
// _zetas = [(mont_r * pow(rou,  br(i,7), q)) % q for i in range(128)]
// zetas = [(int(a)-q) if (a > (q-1)//2) else a for a in _zetas]
const int64_t zetas[PQMX_N] = {
532719943776463611, 1712824210598473412, -1370295354093451833, 124262200749178965, -714223376490581314, -138275710883642869, -1102926763625050, 1767657196266643088, -279197397721253178, 356492044517202178, 1487661701559950979, 1304103426002015710, 1629618216410753009, 1445036730941650804, 1619116696143497796, 105245167723797935, 367723163442095531, -1487071528781225970, 1778506941797319594, -1394647439355336766, 518804722893076255, 141291798759015968, -1463055328214512227, 1403678053431936803, 1629369437266972731, 116935029746966406, 1153404459772611806, 667428405244972246, 1787810736653776682, 347928242243858624, -358756082137868769, -178447815807751549, 998523887703719690, 145591344792672111, -1035093695668908134, -713763388126671157, -1025203798968210061, 918295552878371182, 837927154393850157, -397395448593122216, -749263788221464499, -726376690825334826, 1110731177573170827, -85056957706944670, 1404312931376851606, 400469219639645685, 270310627308674971, 1768608497918996086, 502419640554247474, 931599578710541146, -1180772881415063460, 735294279160574257, 1755306377104268619, -1085009530551386428, 1550895881188125058, -958859323388743482, -50256176979663429, -898615674107773859, 156591195085678140, 1156811988645169641, 716682534159876124, -776091615957528220, 1154435461316943051, 988374064973051516
};


/*************************************************
* Name:        fqmul
*
* Description: Multiplication followed by Montgomery reduction
*
* Arguments:   - int64_t a: first factor
*              - int64_t b: second factor
*
* Returns 16-bit integer congruent to a*b*R^{-1} mod q
**************************************************/
int64_t fqmul(int64_t a, int64_t b) {
  return montgomery_reduce((__int128)a*b);
}

/*************************************************
* Name:        ntt
*
* Description: Inplace number-theoretic transform (NTT) in Rq.
*              input is in standard order, output is in bitreversed order
*
* Arguments:   - int64_t r[4096]: pointer to input/output vector of elements of Zq
**************************************************/
void ntt4(int64_t r[PQMX_N]) {
  unsigned int len, start, j, k;
  int64_t t, zeta;

  k = 1;
  for(len = PQMX_N/2; len >= 4; len >>= 1) {
    for(start = 0; start < PQMX_N; start = j + len) {
      zeta = zetas[k++];
      for(j = start; j < start + len; j++) {
        t = fqmul(zeta, r[j + len]);
        r[j + len] = barrett_reduce(r[j]) - t;
        r[j] = barrett_reduce(r[j]) + t;
      }
    }
  }
}

void ntt8(int64_t r[PQMX_N]) {
  unsigned int len, start, j, k;
  int64_t t, zeta;

  k = 1;
  for(len = PQMX_N/2; len >= 8; len >>= 1) {
    for(start = 0; start < PQMX_N; start = j + len) {
      zeta = zetas[k++];
      for(j = start; j < start + len; j++) {
        t = fqmul(zeta, r[j + len]);
        r[j + len] = barrett_reduce(r[j]) - t;
        r[j] = barrett_reduce(r[j]) + t;
      }
    }
  }
}

void ntt_full(int64_t r[PQMX_N]) {
  unsigned int len, start, j, k;
  int64_t t, zeta;

  k = 1;
  for(len = PQMX_N/2; len > 0; len >>= 1) {
    for(start = 0; start < PQMX_N; start = j + len) {
      zeta = zetas[k++];
      for(j = start; j < start + len; j++) {
        t = fqmul(zeta, r[j + len]);
        r[j + len] = barrett_reduce(r[j]) - t;
        r[j] = barrett_reduce(r[j]) + t;
      }
    }
  }
}

/*************************************************
* Name:        invntt_tomont
*
* Description: Inplace inverse number-theoretic transform in Rq and
*              multiplication by Montgomery factor 2^16.
*              Input is in bitreversed order, output is in standard order
*
* Arguments:   - int16_t r[4096]: pointer to input/output vector of elements of Zq
**************************************************/
void invntt(int64_t r[PQMX_N]) {
  unsigned int start, len, j, k;
  int64_t t, zeta;
  const int64_t f =  PQMX_MONT4; // mont^2/PQMX_L

  k = PQMX_L-1;
  for(len = (PQMX_N/PQMX_L); len < PQMX_N; len <<= 1) {
    for(start = 0; start < PQMX_N; start = j + len) {
      zeta = zetas[k--];
      for(j = start; j < start + len; j++) {
        t = r[j];
        r[j] = barrett_reduce(t + r[j + len]);
        r[j + len] = barrett_reduce(r[j + len] - t);
        r[j + len] = fqmul(zeta, r[j + len]);
      }
    }
  }
  for(j = 0; j < PQMX_N; j++)
    r[j] = fqmul(r[j], f);

}

/*************************************************
* Name:        basemul
*
* Description: Multiplication of polynomials in Zq[X]/(X^4-zeta)
*              used for multiplication of elements in Rq in NTT domain
*
* Arguments:   - int14_t r[4]: pointer to the output polynomial
*              - const int64_t a[4]: pointer to the first factor
*              - const int64_t b[4]: pointer to the second factor
*              - int64_t zeta: integer defining the reduction polynomial
**************************************************/
void basemul4(int64_t r[4], const int64_t a[4], const int64_t b[4], int64_t zeta)
{
  int64_t tmp, rr[4];
  
  rr[0]  = fqmul(a[0], b[0]);
  tmp   = fqmul(a[1], b[3]);
  rr[0] += fqmul(tmp, zeta);
  tmp   = fqmul(a[2], b[2]);
  rr[0]  = barrett_reduce(rr[0]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[3], b[1]);
  rr[0]  = barrett_reduce(rr[0]) + fqmul(tmp, zeta);
  
  rr[1]  = fqmul(a[0], b[1]);
  rr[1]  += fqmul(a[1], b[0]);
  tmp   = fqmul(a[2], b[3]);
  rr[1]  = barrett_reduce(rr[1]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[3], b[2]);
  rr[1]  =  barrett_reduce(rr[1]) + fqmul(tmp, zeta);

  rr[2]  = fqmul(a[0],b[2]);
  rr[2] += fqmul(a[1],b[1]);
  rr[2]  = barrett_reduce(rr[2]) + fqmul(a[2],b[0]);
  tmp   = fqmul(a[3], b[3]);
  rr[2]  = barrett_reduce(rr[2]) + fqmul(tmp, zeta);

  rr[3]  = fqmul(a[0],b[3]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(a[1],b[2]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(a[2],b[1]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(a[3],b[0]);
  
  r[0] = rr[0];
  r[1] = rr[1];
  r[2] = rr[2];
  r[3] = rr[3];
}

void basemul8(int64_t r[8], const int64_t a[8], const int64_t b[8], int64_t zeta)
{
  int64_t tmp, rr[8];

  rr[0]  = fqmul(a[0], b[0]);
  tmp   = fqmul(a[1], b[7]);
  rr[0] += fqmul(tmp, zeta);
  tmp   = fqmul(a[2], b[6]);
  rr[0]  = barrett_reduce(rr[0]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[3], b[5]);
  rr[0]  = barrett_reduce(rr[0]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[4], b[4]);
  rr[0]  = barrett_reduce(rr[0]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[5], b[3]);
  rr[0]  = barrett_reduce(rr[0]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[6], b[2]);
  rr[0]  = barrett_reduce(rr[0]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[7], b[1]);
  rr[0]  = barrett_reduce(rr[0]) + fqmul(tmp, zeta);

  rr[1]  = fqmul(a[0], b[1]);
  rr[1]  += fqmul(a[1], b[0]);
  tmp   = fqmul(a[2], b[7]);
  rr[1]  = barrett_reduce(rr[1]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[3], b[6]);
  rr[1]  =  barrett_reduce(rr[1]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[4], b[5]);
  rr[1]  =  barrett_reduce(rr[1]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[5], b[4]);
  rr[1]  =  barrett_reduce(rr[1]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[6], b[3]);
  rr[1]  =  barrett_reduce(rr[1]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[7], b[2]);
  rr[1]  =  barrett_reduce(rr[1]) + fqmul(tmp, zeta);

  rr[2]  = fqmul(a[0],b[2]);
  rr[2] += fqmul(a[1],b[1]);
  rr[2]  = barrett_reduce(rr[2]) + fqmul(a[2],b[0]);
  tmp   = fqmul(a[3], b[7]);
  rr[2]  = barrett_reduce(rr[2]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[4], b[6]);
  rr[2]  = barrett_reduce(rr[2]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[5], b[5]);
  rr[2]  = barrett_reduce(rr[2]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[6], b[4]);
  rr[2]  = barrett_reduce(rr[2]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[7], b[3]);
  rr[2]  = barrett_reduce(rr[2]) + fqmul(tmp, zeta);

  rr[3]  = fqmul(a[0],b[3]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(a[1],b[2]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(a[2],b[1]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(a[3],b[0]);
  tmp   = fqmul(a[4], b[7]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[5], b[6]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[6], b[5]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[7], b[4]);
  rr[3]  = barrett_reduce(rr[3]) + fqmul(tmp, zeta);

  rr[4]  = fqmul(a[0],b[4]);
  rr[4]  = barrett_reduce(rr[4]) + fqmul(a[1],b[3]);
  rr[4]  = barrett_reduce(rr[4]) + fqmul(a[2],b[2]);
  rr[4]  = barrett_reduce(rr[4]) + fqmul(a[3],b[1]);
  rr[4]  = barrett_reduce(rr[4]) + fqmul(a[4],b[0]);
  tmp   = fqmul(a[5], b[7]);
  rr[4]  = barrett_reduce(rr[4]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[6], b[6]);
  rr[4]  = barrett_reduce(rr[4]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[7], b[5]);
  rr[4]  = barrett_reduce(rr[4]) + fqmul(tmp, zeta);

  rr[5]  = fqmul(a[0],b[5]);
  rr[5]  = barrett_reduce(rr[5]) + fqmul(a[1],b[4]);
  rr[5]  = barrett_reduce(rr[5]) + fqmul(a[2],b[3]);
  rr[5]  = barrett_reduce(rr[5]) + fqmul(a[3],b[2]);
  rr[5]  = barrett_reduce(rr[5]) + fqmul(a[4],b[1]);
  rr[5]  = barrett_reduce(rr[5]) + fqmul(a[5],b[0]);
  tmp   = fqmul(a[6], b[7]);
  rr[5]  = barrett_reduce(rr[5]) + fqmul(tmp, zeta);
  tmp   = fqmul(a[7], b[6]);
  rr[5]  = barrett_reduce(rr[5]) + fqmul(tmp, zeta);

  rr[6]  = fqmul(a[0],b[6]);
  rr[6]  = barrett_reduce(rr[6]) + fqmul(a[1],b[5]);
  rr[6]  = barrett_reduce(rr[6]) + fqmul(a[2],b[4]);
  rr[6]  = barrett_reduce(rr[6]) + fqmul(a[3],b[3]);
  rr[6]  = barrett_reduce(rr[6]) + fqmul(a[4],b[2]);
  rr[6]  = barrett_reduce(rr[6]) + fqmul(a[5],b[1]);
  rr[6]  = barrett_reduce(rr[6]) + fqmul(a[6],b[0]);
  tmp   = fqmul(a[7], b[7]);
  rr[6]  = barrett_reduce(rr[6]) + fqmul(tmp, zeta);

  rr[7]  = fqmul(a[0],b[7]);
  rr[7]  = barrett_reduce(rr[7]) + fqmul(a[1],b[6]);
  rr[7]  = barrett_reduce(rr[7]) + fqmul(a[2],b[5]);
  rr[7]  = barrett_reduce(rr[7]) + fqmul(a[3],b[4]);
  rr[7]  = barrett_reduce(rr[7]) + fqmul(a[4],b[3]);
  rr[7]  = barrett_reduce(rr[7]) + fqmul(a[5],b[2]);
  rr[7]  = barrett_reduce(rr[7]) + fqmul(a[6],b[1]);
  rr[7]  = barrett_reduce(rr[7]) + fqmul(a[7],b[0]);

  r[0] = rr[0];
  r[1] = rr[1];
  r[2] = rr[2];
  r[3] = rr[3];
  r[4] = rr[4];
  r[5] = rr[5];
  r[6] = rr[6];
  r[7] = rr[7];
}

/*************************************************
* Name:        scalar_field_mul
*
* Description: Scalar multiplication of a polynomial in Zq[X]/(X^4-zeta)
*
* Arguments:   - int14_t r[4]: pointer to the output polynomial
*              - const int64_t a: scalar factor
*              - const int64_t b[4]: pointer to the input polynomial factor
**************************************************/
void scalar_field_mul(int64_t r[4], const int64_t a, const int64_t b[4]){
    unsigned int i;
  for(i=0;i<4;i++){
    r[i] = fqmul(a, b[i]);
    r[i] = fqmul(r[i], PQMX_MONT2);
  }
}

/*************************************************
* Name:        field_mul
*
* Description: Field multiplication in Zq.
*
* Arguments:   - int14_t c: pointer to the output element
*              - int64_t a: first field element
*              - int64_t b: second field element
**************************************************/
void field_mul(int64_t *c, int64_t *a, int64_t *b) {
  *c = fqmul(*a, *b);
}