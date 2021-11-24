//MPAMPIS -not included before-
#include "gpu_kernel/blislike.h"

__kernel void blis_like (
    unsigned int k, unsigned int m_alg, unsigned int n_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
  const size_t bx = get_group_id(0);
  const size_t by = get_group_id(1);
  const size_t tx = get_local_id(0);
  const size_t ty = get_local_id(1);

  unsigned int a_offset = bx * BLOCK_SIZE_X;
  unsigned int b_offset = by * BLOCK_SIZE_Y;

  unsigned int ir = bx * BLOCK_SIZE_X / GPU_BLOCK_MR * GPU_BLOCK_MR;
  unsigned int jr = by * BLOCK_SIZE_Y / GPU_BLOCK_NR * GPU_BLOCK_MR;

  unsigned int nr_alg = min((uint) GPU_BLOCK_NR, n_alg - jr);
  unsigned int mr_alg = min((uint) GPU_BLOCK_MR, m_alg - ir);

  unsigned int nr_block = min((uint) BLOCK_SIZE_Y, n_alg - b_offset);
  unsigned int mr_block = min((uint) BLOCK_SIZE_X, m_alg - a_offset);

	__global unsigned int *Ar = &a[ir*k + (a_offset % BLOCK_SIZE_X)];
	__global unsigned int *Br = &b[jr*k + (b_offset % BLOCK_SIZE_Y)];
	__global unsigned int *Cr = &c[a_offset+b_offset*cs_c];

  // can this be shared between threads in a warp??  does it even benefit us?
  __local unsigned int ab[BLOCK_SIZE_Y*BLOCK_SIZE_X];
  __local unsigned int *abij;

  unsigned int bj, ai;
  unsigned int l, i, j;
  for (j = ty; j < nr_block; j += LOCAL_1) {
    abij = &ab[j*mr_block];
    for (i = tx; i < mr_block; i += LOCAL_0) {
      abij[i] = 0;
    }
  }

  // __global unsigned int *Cr_l;
  for (l = 0; l < k; ++l) {
    abij=ab;
    for (j = ty; j < nr_block; j += LOCAL_1) {
      bj = Br[j];      
      abij = &ab[j*mr_block];
      for (i = tx; i < mr_block; i += LOCAL_0) {
        ai = Ar[i];
        abij[i] += popcount(ai & bj);        
      }
    }
    Ar += mr_alg;
    Br += nr_alg;
  }

  __global unsigned int *Cr_l;
  for (j = ty; j < nr_block; j += LOCAL_1) {
    Cr_l = &Cr[j*cs_c];
    abij = &ab[j*mr_block];
    for (i = tx; i < mr_block; i += LOCAL_0) {
      // Cr_l[i*rs_c] += abij[i];
      // since rs_c is assumed to be 1...
      Cr_l[i] += abij[i];
    }
  }
}

__kernel void blis_like_empty (
    unsigned int k, unsigned int mr_alg, unsigned int nr_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
    const size_t by = get_group_id(1);
    const size_t ty = get_local_id(1);  // should be zero?
    const size_t bx = get_group_id(0);
    const size_t tx = get_local_id(0);

    const size_t a_offset = by * GPU_BLOCK_MC + ty;
    const size_t b_offset = bx * GPU_BLOCK_NC + tx;

    // a += a_offset;
    // b += b_offset;
    // is this correct?
    // c += b_offset * cs_c;

    // unsigned int l, i;
    // unsigned int bj, ai;
    // for (l = 0; l < k; ++l) {
    //   bj = ~(*b);
    //   for (i = 0; i < mr_alg; ++i) {
    //     ai = *(a + i);
    //
    //     *(c + i) += popcount(ai & bj);
    //   }
    //   a += mr_alg;
    //   b += nr_alg;
    // }


}

// AxB
__kernel void blis_like8x4 (
    unsigned int k, unsigned int m_alg, unsigned int n_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
  const size_t bx = get_group_id(0);
  const size_t by = get_group_id(1);
  const size_t tx = get_local_id(0);
  const size_t ty = get_local_id(1);

  unsigned int ir = bx * BLOCK_SIZE_X / GPU_BLOCK_MR * GPU_BLOCK_MR;
  unsigned int jr = by * BLOCK_SIZE_Y / GPU_BLOCK_NR * GPU_BLOCK_NR;

  unsigned int mr_alg = min((uint) GPU_BLOCK_MR, m_alg - ir);
  unsigned int nr_alg = min((uint) GPU_BLOCK_NR, n_alg - jr);

  const size_t a_offset = bx * BLOCK_SIZE_X;
  const size_t b_offset = by * BLOCK_SIZE_Y;

  a += ir*k + (a_offset % GPU_BLOCK_MR) + tx;
  b += jr*k + (b_offset % GPU_BLOCK_NR) + ty;
  c += cs_c * (b_offset + ty) + (a_offset + tx);

  unsigned int a0, a1, a2, a3, a4, a5, a6, a7;
  unsigned int b0, b1, b2, b3;

  unsigned int t00, t01, t02, t03;
  unsigned int t10, t11, t12, t13;
  unsigned int t20, t21, t22, t23;
  unsigned int t30, t31, t32, t33;
  unsigned int t40, t41, t42, t43;
  unsigned int t50, t51, t52, t53;
  unsigned int t60, t61, t62, t63;
  unsigned int t70, t71, t72, t73;

  unsigned int c00, c01, c02, c03;
  unsigned int c10, c11, c12, c13;
  unsigned int c20, c21, c22, c23;
  unsigned int c30, c31, c32, c33;
  unsigned int c40, c41, c42, c43;
  unsigned int c50, c51, c52, c53;
  unsigned int c60, c61, c62, c63;
  unsigned int c70, c71, c72, c73;

  c00 = c01 = c02 = c03 = 0;
  c10 = c11 = c12 = c13 = 0;
  c20 = c21 = c22 = c23 = 0;
  c30 = c31 = c32 = c33 = 0;
  c40 = c41 = c42 = c43 = 0;
  c50 = c51 = c52 = c53 = 0;
  c60 = c61 = c62 = c63 = 0;
  c70 = c71 = c72 = c73 = 0;

  unsigned int l;
  for (l = 0; l < k; l++) {

    a0 = *(a + 0*BLOCK_SIZE_X/8);
    a1 = *(a + 1*BLOCK_SIZE_X/8);
    a2 = *(a + 2*BLOCK_SIZE_X/8);
    a3 = *(a + 3*BLOCK_SIZE_X/8);
    a4 = *(a + 4*BLOCK_SIZE_X/8);
    a5 = *(a + 5*BLOCK_SIZE_X/8);
    a6 = *(a + 6*BLOCK_SIZE_X/8);
    a7 = *(a + 7*BLOCK_SIZE_X/8);

    b0 = *(b + 0*BLOCK_SIZE_Y/4);
    b1 = *(b + 1*BLOCK_SIZE_Y/4);
    b2 = *(b + 2*BLOCK_SIZE_Y/4);
    b3 = *(b + 3*BLOCK_SIZE_Y/4);

    a += mr_alg;
    b += nr_alg;

    t00 = a0 & b0;
    t10 = a1 & b0;
    t20 = a2 & b0;
    t30 = a3 & b0;
    t40 = a4 & b0;
    t50 = a5 & b0;
    t60 = a6 & b0;
    t70 = a7 & b0;

    t01 = a0 & b1;
    t11 = a1 & b1;
    t21 = a2 & b1;
    t31 = a3 & b1;
    t41 = a4 & b1;
    t51 = a5 & b1;
    t61 = a6 & b1;
    t71 = a7 & b1;

    t02 = a0 & b2;
    t12 = a1 & b2;
    t22 = a2 & b2;
    t32 = a3 & b2;
    t42 = a4 & b2;
    t52 = a5 & b2;
    t62 = a6 & b2;
    t72 = a7 & b2;

    t03 = a0 & b3;
    t13 = a1 & b3;
    t23 = a2 & b3;
    t33 = a3 & b3;
    t43 = a4 & b3;
    t53 = a5 & b3;
    t63 = a6 & b3;
    t73 = a7 & b3;


    c00 += popcount(t00);
    c10 += popcount(t10);
    c20 += popcount(t20);
    c30 += popcount(t30);
    c40 += popcount(t40);
    c50 += popcount(t50);
    c60 += popcount(t60);
    c70 += popcount(t70);

    c01 += popcount(t01);
    c11 += popcount(t11);
    c21 += popcount(t21);
    c31 += popcount(t31);
    c41 += popcount(t41);
    c51 += popcount(t51);
    c61 += popcount(t61);
    c71 += popcount(t71);

    c02 += popcount(t02);
    c12 += popcount(t12);
    c22 += popcount(t22);
    c32 += popcount(t32);
    c42 += popcount(t42);
    c52 += popcount(t52);
    c62 += popcount(t62);
    c72 += popcount(t72);

    c03 += popcount(t03);
    c13 += popcount(t13);
    c23 += popcount(t23);
    c33 += popcount(t33);
    c43 += popcount(t43);
    c53 += popcount(t53);
    c63 += popcount(t63);
    c73 += popcount(t73);
  }

  *(c + 0*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c00;
  *(c + 1*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c10;
  *(c + 2*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c20;
  *(c + 3*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c30;
  *(c + 4*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c40;
  *(c + 5*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c50;
  *(c + 6*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c60;
  *(c + 7*BLOCK_SIZE_X/8 + cs_c*0*BLOCK_SIZE_Y/4) += c70;

  *(c + 0*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c01;
  *(c + 1*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c11;
  *(c + 2*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c21;
  *(c + 3*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c31;
  *(c + 4*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c41;
  *(c + 5*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c51;
  *(c + 6*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c61;
  *(c + 7*BLOCK_SIZE_X/8 + cs_c*1*BLOCK_SIZE_Y/4) += c71;

  *(c + 0*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c02;
  *(c + 1*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c12;
  *(c + 2*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c22;
  *(c + 3*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c32;
  *(c + 4*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c42;
  *(c + 5*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c52;
  *(c + 6*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c62;
  *(c + 7*BLOCK_SIZE_X/8 + cs_c*2*BLOCK_SIZE_Y/4) += c72;

  *(c + 0*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c03;
  *(c + 1*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c13;
  *(c + 2*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c23;
  *(c + 3*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c33;
  *(c + 4*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c43;
  *(c + 5*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c53;
  *(c + 6*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c63;
  *(c + 7*BLOCK_SIZE_X/8 + cs_c*3*BLOCK_SIZE_Y/4) += c73;
}

// now A is associated with Y, and B with X
__kernel void blis_like4x8 (
    unsigned int k, unsigned int m_alg, unsigned int n_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
  const size_t bx = get_group_id(0);
  const size_t by = get_group_id(1);
  const size_t tx = get_local_id(0);
  const size_t ty = get_local_id(1);

  unsigned int ir = by * BLOCK_SIZE_Y / GPU_BLOCK_MR * GPU_BLOCK_MR;
  unsigned int jr = bx * BLOCK_SIZE_X / GPU_BLOCK_NR * GPU_BLOCK_NR;

  unsigned int mr_alg = min((uint) GPU_BLOCK_MR, m_alg - ir);
  unsigned int nr_alg = min((uint) GPU_BLOCK_NR, n_alg - jr);

  const size_t a_offset = by * BLOCK_SIZE_Y;
  const size_t b_offset = bx * BLOCK_SIZE_X;

  a += ir*k + (a_offset % GPU_BLOCK_MR) + ty;
  b += jr*k + (b_offset % GPU_BLOCK_NR) + tx;
  c += cs_c * (b_offset + tx) + (a_offset + ty);

    unsigned int a0, a1, a2, a3;
    unsigned int b0, b1, b2, b3, b4, b5, b6, b7;

    unsigned int t00, t01, t02, t03, t04, t05, t06, t07;
    unsigned int t10, t11, t12, t13, t14, t15, t16, t17;
    unsigned int t20, t21, t22, t23, t24, t25, t26, t27;
    unsigned int t30, t31, t32, t33, t34, t35, t36, t37;

    unsigned int c00, c01, c02, c03, c04, c05, c06, c07;
    unsigned int c10, c11, c12, c13, c14, c15, c16, c17;
    unsigned int c20, c21, c22, c23, c24, c25, c26, c27;
    unsigned int c30, c31, c32, c33, c34, c35, c36, c37;

    c00 = c01 = c02 = c03 = c04 = c05 = c06 = c07 = 0;
    c10 = c11 = c12 = c13 = c14 = c15 = c16 = c17 = 0;
    c20 = c21 = c22 = c23 = c24 = c25 = c26 = c27 = 0;
    c30 = c31 = c32 = c33 = c34 = c35 = c36 = c37 = 0;

    unsigned int l;
    for (l = 0; l < k; l++) {
      a0 = *(a + 0*BLOCK_SIZE_Y/4);
      a1 = *(a + 1*BLOCK_SIZE_Y/4);
      a2 = *(a + 2*BLOCK_SIZE_Y/4);
      a3 = *(a + 3*BLOCK_SIZE_Y/4);

      b0 = *(b + 0*BLOCK_SIZE_X/8);
      b1 = *(b + 1*BLOCK_SIZE_X/8);
      b2 = *(b + 2*BLOCK_SIZE_X/8);
      b3 = *(b + 3*BLOCK_SIZE_X/8);
      b4 = *(b + 4*BLOCK_SIZE_X/8);
      b5 = *(b + 5*BLOCK_SIZE_X/8);
      b6 = *(b + 6*BLOCK_SIZE_X/8);
      b7 = *(b + 7*BLOCK_SIZE_X/8);

      a += mr_alg;
      b += nr_alg;

      t00 = a0 & b0;
      t01 = a0 & b1;
      t02 = a0 & b2;
      t03 = a0 & b3;
      t04 = a0 & b4;
      t05 = a0 & b5;
      t06 = a0 & b6;
      t07 = a0 & b7;

      t10 = a1 & b0;
      t11 = a1 & b1;
      t12 = a1 & b2;
      t13 = a1 & b3;
      t14 = a1 & b4;
      t15 = a1 & b5;
      t16 = a1 & b6;
      t17 = a1 & b7;

      t20 = a2 & b0;
      t21 = a2 & b1;
      t22 = a2 & b2;
      t23 = a2 & b3;
      t24 = a2 & b4;
      t25 = a2 & b5;
      t26 = a2 & b6;
      t27 = a2 & b7;

      t30 = a3 & b0;
      t31 = a3 & b1;
      t32 = a3 & b2;
      t33 = a3 & b3;
      t34 = a3 & b4;
      t35 = a3 & b5;
      t36 = a3 & b6;
      t37 = a3 & b7;

      c00 += popcount(t00);
      c01 += popcount(t01);
      c02 += popcount(t02);
      c03 += popcount(t03);
      c04 += popcount(t04);
      c05 += popcount(t05);
      c06 += popcount(t06);
      c07 += popcount(t07);

      c10 += popcount(t10);
      c11 += popcount(t11);
      c12 += popcount(t12);
      c13 += popcount(t13);
      c14 += popcount(t14);
      c15 += popcount(t15);
      c16 += popcount(t16);
      c17 += popcount(t17);

      c20 += popcount(t20);
      c21 += popcount(t21);
      c22 += popcount(t22);
      c23 += popcount(t23);
      c24 += popcount(t24);
      c25 += popcount(t25);
      c26 += popcount(t26);
      c27 += popcount(t27);

      c30 += popcount(t30);
      c31 += popcount(t31);
      c32 += popcount(t32);
      c33 += popcount(t33);
      c34 += popcount(t34);
      c35 += popcount(t35);
      c36 += popcount(t36);
      c37 += popcount(t37);
    }

  *(c + cs_c*0*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c00;
  *(c + cs_c*1*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c01;
  *(c + cs_c*2*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c02;
  *(c + cs_c*3*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c03;
  *(c + cs_c*4*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c04;
  *(c + cs_c*5*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c05;
  *(c + cs_c*6*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c06;
  *(c + cs_c*7*BLOCK_SIZE_X/8 + 0*BLOCK_SIZE_Y/4) += c07;

  *(c + cs_c*0*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c10;
  *(c + cs_c*1*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c11;
  *(c + cs_c*2*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c12;
  *(c + cs_c*3*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c13;
  *(c + cs_c*4*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c14;
  *(c + cs_c*5*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c15;
  *(c + cs_c*6*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c16;
  *(c + cs_c*7*BLOCK_SIZE_X/8 + 1*BLOCK_SIZE_Y/4) += c17;

  *(c + cs_c*0*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c20;
  *(c + cs_c*1*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c21;
  *(c + cs_c*2*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c22;
  *(c + cs_c*3*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c23;
  *(c + cs_c*4*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c24;
  *(c + cs_c*5*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c25;
  *(c + cs_c*6*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c26;
  *(c + cs_c*7*BLOCK_SIZE_X/8 + 2*BLOCK_SIZE_Y/4) += c27;

  *(c + cs_c*0*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c30;
  *(c + cs_c*1*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c31;
  *(c + cs_c*2*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c32;
  *(c + cs_c*3*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c33;
  *(c + cs_c*4*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c34;
  *(c + cs_c*5*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c35;
  *(c + cs_c*6*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c36;
  *(c + cs_c*7*BLOCK_SIZE_X/8 + 3*BLOCK_SIZE_Y/4) += c37;

}

// associate X with A, but use a 32x1 local size
__kernel void blis_like4x8v2 (
    unsigned int k, unsigned int m_alg, unsigned int n_alg,
    __global unsigned int *a,
    __global unsigned int *b,
    __global unsigned int *c,
    unsigned int rs_c, unsigned int cs_c
) {
  const size_t bx = get_group_id(0);
  const size_t by = get_group_id(1);
  const size_t tx = get_local_id(0);
  const size_t ty = get_local_id(1);

  unsigned int ir = bx * BLOCK_SIZE_X / GPU_BLOCK_MR * GPU_BLOCK_MR;
  unsigned int jr = by * BLOCK_SIZE_Y / GPU_BLOCK_NR * GPU_BLOCK_NR;

  unsigned int mr_alg = min((uint) GPU_BLOCK_MR, m_alg - ir);
  unsigned int nr_alg = min((uint) GPU_BLOCK_NR, n_alg - jr);

  const size_t a_offset = bx * BLOCK_SIZE_X;
  const size_t b_offset = by * BLOCK_SIZE_Y;

  a += ir*k + (a_offset % GPU_BLOCK_MR) + tx;
  b += jr*k + (b_offset % GPU_BLOCK_NR) + ty;
  c += cs_c * (b_offset + ty) + (a_offset + tx);

    unsigned int a0, a1, a2, a3;
    unsigned int b0, b1, b2, b3, b4, b5, b6, b7;

    unsigned int t00, t01, t02, t03, t04, t05, t06, t07;
    unsigned int t10, t11, t12, t13, t14, t15, t16, t17;
    unsigned int t20, t21, t22, t23, t24, t25, t26, t27;
    unsigned int t30, t31, t32, t33, t34, t35, t36, t37;

    unsigned int c00, c01, c02, c03, c04, c05, c06, c07;
    unsigned int c10, c11, c12, c13, c14, c15, c16, c17;
    unsigned int c20, c21, c22, c23, c24, c25, c26, c27;
    unsigned int c30, c31, c32, c33, c34, c35, c36, c37;

    c00 = c01 = c02 = c03 = c04 = c05 = c06 = c07 = 0;
    c10 = c11 = c12 = c13 = c14 = c15 = c16 = c17 = 0;
    c20 = c21 = c22 = c23 = c24 = c25 = c26 = c27 = 0;
    c30 = c31 = c32 = c33 = c34 = c35 = c36 = c37 = 0;

    unsigned int l;
    for (l = 0; l < k; l++) {
      a0 = *(a + 0*BLOCK_SIZE_X/4);
      a1 = *(a + 1*BLOCK_SIZE_X/4);
      a2 = *(a + 2*BLOCK_SIZE_X/4);
      a3 = *(a + 3*BLOCK_SIZE_X/4);

      b0 = *(b + 0*BLOCK_SIZE_Y/8);
      b1 = *(b + 1*BLOCK_SIZE_Y/8);
      b2 = *(b + 2*BLOCK_SIZE_Y/8);
      b3 = *(b + 3*BLOCK_SIZE_Y/8);
      b4 = *(b + 4*BLOCK_SIZE_Y/8);
      b5 = *(b + 5*BLOCK_SIZE_Y/8);
      b6 = *(b + 6*BLOCK_SIZE_Y/8);
      b7 = *(b + 7*BLOCK_SIZE_Y/8);

      a += mr_alg;
      b += nr_alg;

      t00 = a0 & b0;
      t01 = a0 & b1;
      t02 = a0 & b2;
      t03 = a0 & b3;
      t04 = a0 & b4;
      t05 = a0 & b5;
      t06 = a0 & b6;
      t07 = a0 & b7;

      t10 = a1 & b0;
      t11 = a1 & b1;
      t12 = a1 & b2;
      t13 = a1 & b3;
      t14 = a1 & b4;
      t15 = a1 & b5;
      t16 = a1 & b6;
      t17 = a1 & b7;

      t20 = a2 & b0;
      t21 = a2 & b1;
      t22 = a2 & b2;
      t23 = a2 & b3;
      t24 = a2 & b4;
      t25 = a2 & b5;
      t26 = a2 & b6;
      t27 = a2 & b7;

      t30 = a3 & b0;
      t31 = a3 & b1;
      t32 = a3 & b2;
      t33 = a3 & b3;
      t34 = a3 & b4;
      t35 = a3 & b5;
      t36 = a3 & b6;
      t37 = a3 & b7;

      c00 += popcount(t00);
      c01 += popcount(t01);
      c02 += popcount(t02);
      c03 += popcount(t03);
      c04 += popcount(t04);
      c05 += popcount(t05);
      c06 += popcount(t06);
      c07 += popcount(t07);

      c10 += popcount(t10);
      c11 += popcount(t11);
      c12 += popcount(t12);
      c13 += popcount(t13);
      c14 += popcount(t14);
      c15 += popcount(t15);
      c16 += popcount(t16);
      c17 += popcount(t17);

      c20 += popcount(t20);
      c21 += popcount(t21);
      c22 += popcount(t22);
      c23 += popcount(t23);
      c24 += popcount(t24);
      c25 += popcount(t25);
      c26 += popcount(t26);
      c27 += popcount(t27);

      c30 += popcount(t30);
      c31 += popcount(t31);
      c32 += popcount(t32);
      c33 += popcount(t33);
      c34 += popcount(t34);
      c35 += popcount(t35);
      c36 += popcount(t36);
      c37 += popcount(t37);
    }

  *(c + cs_c*0*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c00;
  *(c + cs_c*1*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c01;
  *(c + cs_c*2*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c02;
  *(c + cs_c*3*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c03;
  *(c + cs_c*4*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c04;
  *(c + cs_c*5*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c05;
  *(c + cs_c*6*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c06;
  *(c + cs_c*7*BLOCK_SIZE_Y/8 + 0*BLOCK_SIZE_X/4) += c07;

  *(c + cs_c*0*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c10;
  *(c + cs_c*1*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c11;
  *(c + cs_c*2*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c12;
  *(c + cs_c*3*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c13;
  *(c + cs_c*4*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c14;
  *(c + cs_c*5*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c15;
  *(c + cs_c*6*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c16;
  *(c + cs_c*7*BLOCK_SIZE_Y/8 + 1*BLOCK_SIZE_X/4) += c17;

  *(c + cs_c*0*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c20;
  *(c + cs_c*1*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c21;
  *(c + cs_c*2*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c22;
  *(c + cs_c*3*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c23;
  *(c + cs_c*4*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c24;
  *(c + cs_c*5*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c25;
  *(c + cs_c*6*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c26;
  *(c + cs_c*7*BLOCK_SIZE_Y/8 + 2*BLOCK_SIZE_X/4) += c27;

  *(c + cs_c*0*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c30;
  *(c + cs_c*1*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c31;
  *(c + cs_c*2*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c32;
  *(c + cs_c*3*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c33;
  *(c + cs_c*4*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c34;
  *(c + cs_c*5*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c35;
  *(c + cs_c*6*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c36;
  *(c + cs_c*7*BLOCK_SIZE_Y/8 + 3*BLOCK_SIZE_X/4) += c37;

}

__kernel void omega (
    __global float *omega, __constant float *ls, __constant float *rs, __constant float *ts,
    __constant int *k, __constant int *m, int inner
) {
  const int i = get_global_id(0);

  const int outer_i = (int)(i/inner);
	const int inner_i = i%inner;

  int vk = k[outer_i];
  int ksel2 = vk * (vk-1) / 2;

  int vm = m[inner_i];
  int msel2 = (vm * (vm-1)) / 2;

  float numerator = (ls[outer_i] + rs[inner_i]) / (ksel2 + msel2);

  float denominator = (ts[i] - ls[outer_i] - rs[inner_i]) / (k[outer_i]*m[inner_i]) + 0.00001f;

  omega[i] =  numerator / denominator;

  // omega[get_global_id(0)] = 33.21;
  // omega[get_global_id(0)] = ts[get_global_id(0)];
  // omega[i] = ts[i];  // Using i is faster than twice "get_global_id(0)"
}

__kernel void omega2 (
    __global float *omega, __constant float *ls, __constant float *rs, __constant float *ts,
    __constant int *k, __constant int *m
) {
  const int i = get_global_id(0);

  int ksel2 = k[i] * (k[i]-1) / 2;

  int msel2 = (m[i] * (m[i]-1)) / 2;

  float numerator = (ls[i] + rs[i]) / (ksel2 + msel2);

  float denominator = (ts[i] - ls[i] - rs[i]) / (k[i]*m[i]) + 0.00001f;

  omega[i] =  numerator / denominator;

  // Make use of __local omega array and do a barrier sync and then a for loop for finding max
  // Also see AMD guide for filling kernel or hiding mem overhead with alu operations
  // Use int4 and float4??
}

__kernel void omega3 (
    __global float *omega, __constant float *lr, __constant float *ts, __constant int *km, int inner
) {
  const unsigned int i = get_global_id(0);

  unsigned int outer_i = i / inner + inner;
	unsigned int inner_i = i % inner;

  int vk = km[outer_i];
  int ksel2 = (vk * (vk-1)) / 2;

  int vm = km[inner_i];
  int msel2 = (vm * (vm-1)) / 2;

  float n = (lr[outer_i] + lr[inner_i]) / (ksel2 + msel2);

  float d = (ts[i] - lr[outer_i] - lr[inner_i]) / (vk*vm) + 0.00001f;

  omega[i] = n / d;
}

// __kernel void omega4 (
//     __global float *omega_global, __local float *omega_local, __global unsigned int *index_global, 
//     __local unsigned int *index_local, __constant float *lrkm, __constant float *ts, int inner, 
//     unsigned int iter
// ) {
//   unsigned int ig = get_global_id(0);
//   unsigned int il = get_local_id(0);
//   unsigned int size_wg = get_local_size(0);
//   unsigned int size_cu = get_num_groups(0);
//   unsigned int iwg = get_group_id(0);
//   unsigned int ic = ig * iter;

//   float numerator, denominator, maxW=0.0f, tmpW;
//   unsigned int maxI, l, id, io, ii;
//   int vk, ksel, vm, msel;
  
//   for(l = 0; l < iter; l++){    // Read "4 bytes/cycle" more reads higher speeds
//     id = ic + l;
//     io = (int)(id / inner) * 2 + (2 * inner);//((int)(id / inner) << 1) + (inner << 1);
//     ii = (id % inner) * 2;//(id % inner) << 1;

//     vk = (int)lrkm[io+1];
//     ksel = (vk * (vk-1)) / 2;//(vk * (vk-1)) >> 1;

//     vm = (int)lrkm[ii+1];
//     msel = (vm * (vm-1)) / 2;//(vm * (vm-1)) >> 1;

//     numerator = (lrkm[io] + lrkm[ii]) / (ksel + msel);

//     denominator = (ts[id] - lrkm[io] - lrkm[ii]) / (vk*vm) + 0.00001f;

//     tmpW = numerator / denominator;

//     if(tmpW > maxW){
//       maxW = tmpW;
//       maxI = id;
//     }
//   }
//   omega_local[il] = maxW;
//   index_local[il] = maxI;
//   // Wait for all local threads to finish 
//   barrier(CLK_LOCAL_MEM_FENCE);
//   if(il == 0){  // Every 1st thread of a work group // Can use multiple threads with tmpW2, maxW2
//     for(l = 1; l < size_wg; l++){   // For every local memory item (work-group item)
//       if(omega_local[l] > omega_local[0]){  // number of work groups remaining maxW values
//         omega_local[0] = omega_local[l];
//         index_local[0] = index_local[l];
//       }
//     }
//     omega_global[iwg] = omega_local[0];
//     index_global[iwg] = index_local[0];
//   }
//   barrier(CLK_GLOBAL_MEM_FENCE);
//   if(ig == 0){  // Very first (global) thread
//     for(l = 1; l < size_cu; l++){   // For every global memory item (work-group id)
//       if(omega_global[l] > omega_global[0]){
//         omega_global[0] = omega_global[l];
//         index_global[0] = index_global[l];
//       }
//     }
//   }
// }

// __kernel void omega5 (
//     __global float *omega_global, __local float *omega_local, __global unsigned int *index_global, 
//     __local unsigned int *index_local, __constant float *ls, __constant float *rs, __constant float *ts,
//     __constant int *k, __constant int *m, unsigned int iter
// ) {
//   unsigned int ig = get_global_id(0);
//   unsigned int il = get_local_id(0);
//   unsigned int size_wg = get_local_size(0);
//   unsigned int size_cu = get_num_groups(0);
//   unsigned int iwg = get_group_id(0);
//   unsigned int ic = ig * iter;

//   float numerator, denominator, maxW=0.0f, tmpW;
//   unsigned int maxI, l, id, io, ii;
//   int ksel, msel;
  
//   for(l = 0; l < iter; l++){    // Read "4 bytes/cycle" more reads higher speeds? Bank?
//     id = ic + l;

//     ksel = k[id] * (k[id]-1) / 2;

//     msel = (m[id] * (m[id]-1)) / 2;

//     numerator = (ls[id] + rs[id]) / (ksel + msel);

//     denominator = (ts[id] - ls[id] - rs[id]) / (k[id]*m[id]) + 0.00001f;

//     tmpW = numerator / denominator;

//     if(tmpW > maxW){
//       maxW = tmpW;
//       maxI = id;
//     }
//   }
//   omega_local[il] = maxW;
//   index_local[il] = maxI;
//   // Wait for all local threads to finish 
//   barrier(CLK_LOCAL_MEM_FENCE);
//   if(il == 0){  // Every 1st thread of a work group // Can use multiple threads with tmpW2, maxW2
//     for(l = 1; l < size_wg; l++){   // For every local memory item (work-group item)
//       if(omega_local[l] > omega_local[0]){  // number of work groups remaining maxW values
//         omega_local[0] = omega_local[l];
//         index_local[0] = index_local[l];
//       }
//     }
//     omega_global[iwg] = omega_local[0];
//     index_global[iwg] = index_local[0];
//   }
//   barrier(CLK_GLOBAL_MEM_FENCE);
//   if(ig == 0){  // Very first (global) thread
//     for(l = 1; l < size_cu; l++){   // For every global memory item (work-group id)
//       if(omega_global[l] > omega_global[0]){
//         omega_global[0] = omega_global[l];
//         index_global[0] = index_global[l];
//       }
//     }
//   }
// }

__kernel void omega6 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __constant float *ts, __constant int *km, int inner, unsigned int iter
) {
  unsigned int ig = get_global_id(0);
  unsigned int ic = ig * iter;

  float numerator, denominator, maxW=0.0f, tmpW;
  unsigned int maxI, i, id, io, ii;
  int vk, ksel, vm, msel;
  
  for(i = 0; i < iter; i++){
    id = ic + i;
    io = id / inner + inner;
    ii = id % inner;

    vk = km[io];
    ksel = (vk * (vk-1)) / 2;//(vk * (vk-1)) >> 1;

    vm = km[ii];
    msel = (vm * (vm-1)) / 2;//(vm * (vm-1)) >> 1;

    numerator = (lr[io] + lr[ii]) / (ksel + msel);

    denominator = (ts[id] - lr[io] - lr[ii]) / (vk*vm) + 0.00001f;

    tmpW = numerator / denominator;

    if(tmpW > maxW){
      maxW = tmpW;
      maxI = id;
    }
  }
  omega_global[ig] = maxW;
  // omega_global[ig+1] = (float)maxI;
  index_global[ig] = maxI;
}

__kernel void omega7 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __constant float *ts, __constant int *km, int inner, unsigned int iter
) {
  unsigned int ig = get_global_id(0);
  unsigned int ic = ig * iter;

  const float den_off = 0.00001f;
  float maxW=0.0f;
  unsigned int maxI, i, id, lf = 4;

  float l1, l2, l3, l4;
  float r1, r2, r3, r4;
  float t1, t2, t3, t4;
  float n1, n2, n3, n4;
  float d1, d2, d3, d4;
  float tmp1, tmp2, tmp3, tmp4;

  int k1, k2, k3, k4;
  int m1, m2, m3, m4;
  int ks1, ks2, ks3, ks4;
  int ms1, ms2, ms3, ms4;

  unsigned int io1, io2, io3, io4;
  unsigned int ii1, ii2, ii3, ii4;

  // io = (ic / inner) + inner;
  // ii = ic % inner;
  // rest = inner - ii - 4;
  // bool done = true;
  // if(rest < iter)
  //   done = false;

  for(i = 0; i+lf-1 < iter; i+=lf){
    id = ic + i;
    io1 = (id / inner) + inner;
    io2 = ((id + 1) / inner) + inner;
    io3 = ((id + 2) / inner) + inner;
    io4 = ((id + 3) / inner) + inner;

    ii1 = id % inner;
    ii2 = (id + 1) % inner;
    ii3 = (id + 2) % inner;
    ii4 = (id + 3) % inner;

    l1 = lr[io1];
    l2 = lr[io2];
    l3 = lr[io3];
    l4 = lr[io4];

    k1 = km[io1];
    k2 = km[io2];
    k3 = km[io3];
    k4 = km[io4];

    r1 = lr[ii1];
    r2 = lr[ii2];
    r3 = lr[ii3];
    r4 = lr[ii4];
    
    m1 = km[ii1];
    m2 = km[ii2];
    m3 = km[ii3];
    m4 = km[ii4];
    
    t1 = ts[id];
    t2 = ts[id + 1];
    t3 = ts[id + 2];
    t4 = ts[id + 3];

    /*if(i > rest){
      int index1 = (ic    ) / inner + inner;
      int index2 = (ic + 1) / inner + inner;
      int index3 = (ic + 2) / inner + inner;
      int index4 = (ic + 3) / inner + inner;
      l1 = lr[index1];
      l2 = lr[index2];
      l3 = lr[index3];
      l4 = lr[index4];

      k1 = km[index1];
      k2 = km[index2];
      k3 = km[index3];
      k4 = km[index4];
      
      int index11 = (ic    ) % inner;
      int index21 = (ic + 1) % inner;
      int index31 = (ic + 2) % inner;
      int index41 = (ic + 3) % inner;

      r1 = lr[index11];
      r2 = lr[index21];
      r3 = lr[index31];
      r4 = lr[index41];

      m1 = km[index11];
      m2 = km[index21];
      m3 = km[index31];
      m4 = km[index41];

      io++;
      ii = i - rest;
    }
    ii++;
    ic++;
    */
    ks1 = (k1 * (k1-1)) / 2;
    ks2 = (k2 * (k2-1)) / 2;
    ks3 = (k3 * (k3-1)) / 2;
    ks4 = (k4 * (k4-1)) / 2;

    ms1 = (m1 * (m1-1)) / 2;
    ms2 = (m2 * (m2-1)) / 2;
    ms3 = (m3 * (m3-1)) / 2;
    ms4 = (m4 * (m4-1)) / 2;

    n1 = (l1 + r1) / (ks1 + ms1);
    n2 = (l2 + r2) / (ks2 + ms2);
    n3 = (l3 + r3) / (ks3 + ms3);
    n4 = (l4 + r4) / (ks4 + ms4);

    d1 = (t1 - l1 - r1) / (k1 * m1) + den_off;
    d2 = (t2 - l2 - r2) / (k2 * m2) + den_off;
    d3 = (t3 - l3 - r3) / (k3 * m3) + den_off;
    d4 = (t4 - l4 - r4) / (k4 * m4) + den_off;

    tmp1 = n1 / d1;
    tmp2 = n2 / d2;
    tmp3 = n3 / d3;
    tmp4 = n4 / d4;

    if(tmp1 > maxW){
      maxW = tmp1;
      maxI = id;
    }
    if(tmp2 > maxW){
      maxW = tmp2;
      maxI = id + 1;
    }
    if(tmp3 > maxW){
      maxW = tmp3;
      maxI = id + 2;
    }
    if(tmp4 > maxW){
      maxW = tmp4;
      maxI = id + 3;
    }
  }
  for(/*emp*/;i<iter;i++){
    id = ic + i;
    io1 = (id / inner) + inner;
    ii1 = id % inner;

    l1 = lr[io1];
    k1 = km[io1];
    r1 = lr[ii1];
    m1 = km[ii1];
    t1 = ts[id];

    ks1 = (k1 * (k1-1)) / 2;

    ms1 = (m1 * (m1-1)) / 2;

    n1 = (l1 + r1) / (ks1 + ms1);

    d1 = (t1 - l1 - r1) / (k1 * m1) + den_off;

    tmp1 = n1 / d1;
    
    if(tmp1 > maxW){
      maxW = tmp1;
      maxI = id;
    }
  }
  omega_global[ig] = maxW;
  index_global[ig] = maxI;
}

__kernel void omega8 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __constant float *ts, __constant int *km, int inner
) {
  unsigned int ig = get_global_id(0);
  unsigned int io = ig + inner;
  unsigned int ic = ig * inner;

  const float den_off = 0.00001f;
  unsigned int maxI, i, id;

  float l, r, t, n, d, tmpW, maxW = 0.0f;
  int k, m, ks, ms;

  l = lr[io];
  k = km[io];
  ks = (k * (k-1)) / 2;

  for(i = 0; i < inner; i++){
    id = ic + i;
    r = lr[i];
    m = km[i];
    
    t = ts[id];
    // t = ts[ic];

    ms = (m * (m-1)) / 2;
    n = (l + r) / (ks + ms);
    d = (t - l - r) / (k * m) + den_off;
    tmpW = n / d;

    if(tmpW > maxW){
      maxW = tmpW;
      maxI = id;
      // maxI = ic;
    }
    // ic++;
  }
  omega_global[ig] = maxW;
  index_global[ig] = maxI;
}

__kernel void omega9 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __constant float *ts, __constant int *km, int inner
) {
  unsigned int ig = get_global_id(0);
  unsigned int io = ig + inner;
  unsigned int ic = ig * inner;
  // unsigned int ie = ic + inner;

  const float den_off = 0.00001f;
  float maxW=0.0f;
  unsigned int maxI, i, id, lf = 4;

  float l1;
  float r1, r2, r3, r4;
  float t1, t2, t3, t4;
  float n1, n2, n3, n4;
  float d1, d2, d3, d4;
  float tmp1, tmp2, tmp3, tmp4;

  int k1;
  int m1, m2, m3, m4;
  int ks1;
  int ms1, ms2, ms3, ms4;

  unsigned int ii = 0;

  l1 = lr[io];
  k1 = km[io];

  ks1 = (k1 * (k1-1)) / 2;

  for(i = 0; i+lf-1 < inner; i+=lf){
    id = ic + i;
    r1 = lr[ii];
    m1 = km[ii++];
    r2 = lr[ii];
    m2 = km[ii++];
    r3 = lr[ii];
    m3 = km[ii++];
    r4 = lr[ii];
    m4 = km[ii++];
    
    t1 = ts[id];
    t2 = ts[id + 1];
    t3 = ts[id + 2];
    t4 = ts[id + 3];

    ms1 = (m1 * (m1-1)) / 2;
    ms2 = (m2 * (m2-1)) / 2;
    ms3 = (m3 * (m3-1)) / 2;
    ms4 = (m4 * (m4-1)) / 2;

    n1 = (l1 + r1) / (ks1 + ms1);
    n2 = (l1 + r2) / (ks1 + ms2);
    n3 = (l1 + r3) / (ks1 + ms3);
    n4 = (l1 + r4) / (ks1 + ms4);

    d1 = (t1 - l1 - r1) / (k1 * m1) + den_off;
    d2 = (t2 - l1 - r2) / (k1 * m2) + den_off;
    d3 = (t3 - l1 - r3) / (k1 * m3) + den_off;
    d4 = (t4 - l1 - r4) / (k1 * m4) + den_off;

    tmp1 = n1 / d1;
    tmp2 = n2 / d2;
    tmp3 = n3 / d3;
    tmp4 = n4 / d4;

    if(tmp1 > maxW){
      maxW = tmp1;
      maxI = id;
    }
    if(tmp2 > maxW){
      maxW = tmp2;
      maxI = id + 1;
    }
    if(tmp3 > maxW){
      maxW = tmp3;
      maxI = id + 2;
    }
    if(tmp4 > maxW){
      maxW = tmp4;
      maxI = id + 3;
    }
  }
  for(/*emp*/;i<inner;i++){
    id = ic + i;
    r1 = lr[ii];
    m1 = km[ii++];
    t1 = ts[id];

    ms1 = (m1 * (m1-1)) / 2;

    n1 = (l1 + r1) / (ks1 + ms1);

    d1 = (t1 - l1 - r1) / (k1 * m1) + den_off;

    tmp1 = n1 / d1;
    
    if(tmp1 > maxW){
      maxW = tmp1;
      maxI = id;
    }
  }
  omega_global[ig] = maxW;
  index_global[ig] = maxI;
}

__kernel void omega10 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __constant float *ts, __constant int *km, int inner
) {
  unsigned int ig = get_global_id(0);
  unsigned int io = ig + inner;
  unsigned int ic = ig * inner;

  const float den_off = 0.00001f;
  unsigned int maxI, i, id;

  float l, r, t, n, d, tmpW, maxW = 0.0f;
  int k, m, ks, ms;

  l = lr[io];
  k = km[io];
  ks = (k * (k-1)) / 2;

  for(i = 0; i < inner; i++){
    id = ic + i;
    r = lr[i];
    m = km[i];
    
    t = ts[id];
    // t = ts[ic];

    ms = (m * (m-1)) / 2;
    n = (l + r) / (ks + ms);
    d = (t - l - r) / (k * m) + den_off;
    tmpW = n / d;

    if(tmpW > maxW){
      maxW = tmpW;
      maxI = id;
      // maxI = ic;
    }
    // ic++;
  }
  omega_global[ig] = maxW;
  index_global[ig] = maxI;
}

__kernel void omegatest (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __local float *r_local, __constant float *ts, __constant int *km, __local int *m_local, int inner
) {
  unsigned int ig = get_global_id(0);
  unsigned int gs = get_global_size(0);
  unsigned int il = get_local_id(0);
  unsigned int wg = get_group_id(0);
  unsigned int ls = get_local_size(0);
  unsigned int ic = ig * inner;

  unsigned int i, id, im = ls * inner * wg, ip = 0;
  float t, tmp = 0.0f;

  for(i = 0; i < inner; i++){
  // for(i = ic; i < ic+inner; i++){
  //   t = ts[i];

    // t = ts[im + il];
    // im+=ls;

    t = ts[ip + ig];
    ip+=gs;

    tmp += t;
  }
  // tmp += ts[ig];
  omega_global[ig] = tmp;
}

__kernel void omega12 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __local float *r_local, __constant float *ts, __constant int *km, __local int *m_local, int inner
) {
  unsigned int ig = get_global_id(0);
  unsigned int gs = get_global_size(0);
  unsigned int il = get_local_id(0);
  unsigned int ws = get_local_size(0);
  unsigned int io = ig + inner;

  const float den_off = 0.00001f;
  unsigned int maxI;
  unsigned int o, p, q;

  float l, r, t, n, d, tmpW, maxW = 0.0f;
  int k, m, ks, ms;

  l = lr[io];
  k = km[io];
  ks = (k * (k-1)) / 2;

  q = inner / ws;

  for(o = 0; o < q; o++){
    r_local[il] = lr[ws * o + il];
    m_local[il] = km[ws * o + il];
    // Ensure writes have completed:
    // barrier(CLK_LOCAL_MEM_FENCE);   // Is stated in optimization guide since work group size is > 64
    for(p = 0; p < ws; p++){
      r = r_local[p];
      m = m_local[p];
      t = ts[(ws * o + p) * gs + ig];

      ms = (m * (m-1)) / 2;
      n = (l + r) / (ks + ms);
      d = (t - l - r) / (k * m) + den_off;
      tmpW = n / d;

      if(tmpW > maxW){
        maxW = tmpW;
        maxI = ws * n + p;
      }
    }
  }

  omega_global[ig] = maxW;
  index_global[ig] = maxI;
}

__kernel void omega13 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __constant float *ts, __constant int *km, int mult, int iter, int inner
) {
  unsigned int ig = get_global_id(0);
  unsigned int gs = get_global_size(0);
  unsigned int outer = gs / mult;
  unsigned int io = ig / mult;
  unsigned int is = (ig % mult * iter) + outer;
  unsigned int ic = io * inner - outer;

  const float den_off = 0.00001f;
  unsigned int maxI, i, ip = 0;

  float l, r, t, n, d, tmpW, maxW = 0.0f;
  int k, m, ks, ms;

  l = lr[io];
  k = km[io];
  // l = lr[ig];
  // k = km[ig];
  ks = (k * (k-1)) / 2;
  #pragma unroll 8
  for(i = is; i < iter+is; i++){
    r = lr[i];
    m = km[i];

    t = ts[ic + i];
    // iss++;
    
    // t = ts[ip + ig];
    // ip += gs;

    ms = (m * (m-1)) / 2;
    n = (l + r) / (ks + ms);
    d = (t - l - r) / (k * m) + den_off;
    tmpW = n / d;

    if(tmpW > maxW){
      maxW = tmpW;
      maxI = i;
    }
  }
  omega_global[ig] = maxW;
  index_global[ig] = maxI + ic;
}

__kernel void omega14 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *ls, 
    __constant float *rs, __constant float *ts, __constant int *kss, __constant int *mss, int mult,
    int iter, int inner
) {
  unsigned int ig = get_global_id(0);
  unsigned int gs = get_global_size(0);
  unsigned int outer = gs / mult;
  unsigned int io = ig / mult;
  unsigned int is = (ig % mult * iter);
  unsigned int ic = io * inner;

  const float den_off = 0.00001f;
  unsigned int maxI, i, ip = 0;
  // unsigned int maxI, i, lf = 4, ip = ig * lf;
  // unsigned int maxI1, maxI2, maxI3, maxI4;

  // float tmpW, maxW = 0.0f;
  float l, r, t, n, d, tmpW, maxW = 0.0f;
  int k, m, ks, ms;

  // float l1;
  // float r1, r2, r3, r4, r5, r6, r7, r8;
  // float t1, t2, t3, t4, t5, t6, t7, t8;
  // float n1, n2, n3, n4, n5, n6, n7, n8;
  // float d1, d2, d3, d4, d5, d6, d7, d8;
  // float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  // float maxW1=0.0f, maxW2=0.0f, maxW3=0.0f, maxW4=0.0f;

  // int k1;
  // int m1, m2, m3, m4, m5, m6, m7, m8;
  // int ks1;
  // int ms1, ms2, ms3, ms4, ms5, ms6, ms7, ms8;

  l = ls[io];
  k = kss[io];
  // l = ls[ig];
  // k = kss[ig];
  ks = (k * (k-1)) / 2;

  #pragma unroll 4
  for(i = is; i < iter + is; i++){
  // for(i = 0; i < iter; i++){
    r = rs[ip + ig];
    m = mss[ip + ig];
    // r = rs[i];
    // m = mss[i];
    t = ts[ip + ig];
    ip += gs;

    ms = (m * (m-1)) / 2;
    n = (l + r) / (ks + ms);
    d = (t - l - r) / (k * m) + den_off;
    tmpW = n / d;

    if(tmpW > maxW){
      maxW = tmpW;
      maxI = i;
    }
  }

  // l1 = ls[io];
  // k1 = kss[io];
  // ks1 = (k1 * (k1-1)) / 2;

  // for(i = 0; i+lf-1 < iter; i+=lf){
  //   r1 = rs[ip + ig    ];
  //   r2 = rs[ip + ig + 1];
  //   r3 = rs[ip + ig + 2];
  //   r4 = rs[ip + ig + 3];
  //   // r5 = rs[ip + ig + 4];
  //   // r6 = rs[ip + ig + 5];
  //   // r7 = rs[ip + ig + 6];
  //   // r8 = rs[ip + ig + 7];
  //   m1 = mss[ip + ig    ];
  //   m2 = mss[ip + ig + 1];
  //   m3 = mss[ip + ig + 2];
  //   m4 = mss[ip + ig + 3];
  //   // m5 = mss[ip + ig + 4];
  //   // m6 = mss[ip + ig + 5];
  //   // m7 = mss[ip + ig + 6];
  //   // m8 = mss[ip + ig + 7];
  //   t1 = ts[ip + ig    ];
  //   t2 = ts[ip + ig + 1];
  //   t3 = ts[ip + ig + 2];
  //   t4 = ts[ip + ig + 3];
  //   // t5 = ts[ip + ig + 4];
  //   // t6 = ts[ip + ig + 5];
  //   // t7 = ts[ip + ig + 6];
  //   // t8 = ts[ip + ig + 7];
  //   ip += lf * gs;

  //   ms1 = (m1 * (m1-1)) * 0.5;
  //   ms2 = (m2 * (m2-1)) * 0.5;
  //   ms3 = (m3 * (m3-1)) * 0.5;
  //   ms4 = (m4 * (m4-1)) * 0.5;
  //   // ms5 = (m5 * (m5-1)) / 2;
  //   // ms6 = (m6 * (m6-1)) / 2;
  //   // ms7 = (m7 * (m7-1)) / 2;
  //   // ms8 = (m8 * (m8-1)) / 2;

  //   n1 = (l1 + r1) / (ks1 + ms1);
  //   n2 = (l1 + r2) / (ks1 + ms2);
  //   n3 = (l1 + r3) / (ks1 + ms3);
  //   n4 = (l1 + r4) / (ks1 + ms4);
  //   // n5 = (l1 + r5) / (ks1 + ms5);
  //   // n6 = (l1 + r6) / (ks1 + ms6);
  //   // n7 = (l1 + r7) / (ks1 + ms7);
  //   // n8 = (l1 + r8) / (ks1 + ms8);

  //   d1 = (t1 - l1 - r1) / (k1 * m1) + den_off;
  //   d2 = (t2 - l1 - r2) / (k1 * m2) + den_off;
  //   d3 = (t3 - l1 - r3) / (k1 * m3) + den_off;
  //   d4 = (t4 - l1 - r4) / (k1 * m4) + den_off;
  //   // d5 = (t5 - l1 - r5) / (k1 * m5) + den_off;
  //   // d6 = (t6 - l1 - r6) / (k1 * m6) + den_off;
  //   // d7 = (t7 - l1 - r7) / (k1 * m7) + den_off;
  //   // d8 = (t8 - l1 - r8) / (k1 * m8) + den_off;

  //   tmp1 = n1 / d1;
  //   tmp2 = n2 / d2;
  //   tmp3 = n3 / d3;
  //   tmp4 = n4 / d4;
  //   // tmp5 = n5 / d5;
  //   // tmp6 = n6 / d6;
  //   // tmp7 = n7 / d7;
  //   // tmp8 = n8 / d8;

  //   if(tmp1 > maxW){
  //     maxW = tmp1;
  //     maxI = i;
  //   }
  //   if(tmp2 > maxW){
  //     maxW = tmp2;
  //     maxI = i + 1;
  //   }
  //   if(tmp3 > maxW){
  //     maxW = tmp3;
  //     maxI = i + 2;
  //   }
  //   if(tmp4 > maxW){
  //     maxW = tmp4;
  //     maxI = i + 3;
  //   }
  //   // if(tmp5 > maxW){
  //   //   maxW = tmp5;
  //   //   maxI = ic + 4;
  //   // }
  //   // if(tmp6 > maxW){
  //   //   maxW = tmp6;
  //   //   maxI = ic + 5;
  //   // }
  //   // if(tmp7 > maxW){
  //   //   maxW = tmp7;
  //   //   maxI = ic + 6;
  //   // }
  //   // if(tmp8 > maxW){
  //   //   maxW = tmp8;
  //   //   maxI = ic + 7;
  //   // }
  // }
  // // maxW = maxW1;
  // // if(maxW2 > maxW)
  // //   maxW = maxW2;
  // // if(maxW3 > maxW)
  // //   maxW = maxW3;
  // // if(maxW4 > maxW)
  // //   maxW = maxW4;

  // for(/*emp*/;i<iter;i++){
  //   r1 = rs[ip + ig];
  //   m1 = mss[ip + ig];
  //   t1 = ts[ip + ig];
  //   ip += gs;

  //   ms1 = (m1 * (m1-1)) / 2;

  //   n1 = (l1 + r1) / (ks1 + ms1);

  //   d1 = (t1 - l1 - r1) / (k1 * m1) + den_off;

  //   tmp1 = n1 / d1;
    
  //   if(tmp1 > maxW){
  //     maxW = tmp1;
  //     maxI = i;
  //   }
  // }
  omega_global[ig] = maxW;
  index_global[ig] = maxI + ic;
}

// __kernel void omega15 (
//     __global float *omega_global, __global unsigned int *index_global, __constant float *ls, 
//     __constant float *rs, __constant float *ts, __constant int *kss, __constant int *mss, int gr_load,
//     int it_load
// ) {
//   unsigned int wg = get_group_id(0);
//   unsigned int il = get_local_id(0);
//   unsigned int gs = get_local_size(0);
//   unsigned int ig = get_global_id(0);
//   unsigned int gr = get_num_groups(0);
//   // unsigned int ii = wg * gs * (gr_load * it_load) + il;
//   unsigned int ii = ig; //wg * gs + il;

//   const float den_off = 0.00001f;
//   unsigned int maxI, i, j;
//   // unsigned int maxI, i, j, lf = 8;

//   // float tmpW, maxW = 0.0f;
//   float l, r, t, n, d, tmpW, maxW = 0.0f;
//   int k, m, ks, ms;

//   // float l1;
//   // float r1, r2, r3, r4, r5, r6, r7, r8;
//   // float t1, t2, t3, t4, t5, t6, t7, t8;
//   // float n1, n2, n3, n4, n5, n6, n7, n8;
//   // float d1, d2, d3, d4, d5, d6, d7, d8;
//   // float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

//   // int k1;
//   // int m1, m2, m3, m4, m5, m6, m7, m8;
//   // int ks1;
//   // int ms1, ms2, ms3, ms4, ms5, ms6, ms7, ms8;

//   // #pragma unroll
//   for(i = 0; i < gr_load; i++){
//     l = ls[wg];         // wg = wgID + (i * #wg);
//     k = kss[wg];
//     wg += gr;
//     ks = (k * (k-1)) / 2;
//     // #pragma unroll 8
//     for(j = 0; j < it_load; j++){
//       r = rs[ii];       // ii = (wgID * wgSIZE) + (j + (i * it_load) * #wg * wgSIZE) + localID;
//       m = mss[ii];
//       t = ts[ii];

//       ms = (m * (m-1)) / 2;
//       n = (l + r) / (ks + ms);
//       d = (t - l - r) / (k * m) + den_off;
//       tmpW = n / d;

//       if(tmpW > maxW){
//         maxW = tmpW;
//         maxI = ii;
//       }
//       // ii += gs;
//       ii += gs * gr;
//     }
//   // for(i = 0; i < gr_load; i++){
//   //   l1 = ls[wg];         // wg = wgID + (i * #wg);
//   //   k1 = kss[wg];
//   //   wg += gr;
//   //   ks1 = (k1 * (k1-1)) / 2;
//   //   for(j = 0; j+lf-1 < it_load; j+=lf){
//   //     r1 = rs[ii    ];
//   //     r2 = rs[ii + 1];
//   //     r3 = rs[ii + 2];
//   //     r4 = rs[ii + 3];
//   //     r5 = rs[ii + 4];
//   //     r6 = rs[ii + 5];
//   //     r7 = rs[ii + 6];
//   //     r8 = rs[ii + 7];
//   //     m1 = mss[ii    ];
//   //     m2 = mss[ii + 1];
//   //     m3 = mss[ii + 2];
//   //     m4 = mss[ii + 3];
//   //     m5 = mss[ii + 4];
//   //     m6 = mss[ii + 5];
//   //     m7 = mss[ii + 6];
//   //     m8 = mss[ii + 7];
//   //     t1 = ts[ii    ];
//   //     t2 = ts[ii + 1];
//   //     t3 = ts[ii + 2];
//   //     t4 = ts[ii + 3];
//   //     t5 = ts[ii + 4];
//   //     t6 = ts[ii + 5];
//   //     t7 = ts[ii + 6];
//   //     t8 = ts[ii + 7];
//   //     ii += lf * gs * gr;

//   //     ms1 = (m1 * (m1-1)) / 2;
//   //     ms2 = (m2 * (m2-1)) / 2;
//   //     ms3 = (m3 * (m3-1)) / 2;
//   //     ms4 = (m4 * (m4-1)) / 2;
//   //     ms5 = (m5 * (m5-1)) / 2;
//   //     ms6 = (m6 * (m6-1)) / 2;
//   //     ms7 = (m7 * (m7-1)) / 2;
//   //     ms8 = (m8 * (m8-1)) / 2;

//   //     n1 = (l1 + r1) / (ks1 + ms1);
//   //     n2 = (l1 + r2) / (ks1 + ms2);
//   //     n3 = (l1 + r3) / (ks1 + ms3);
//   //     n4 = (l1 + r4) / (ks1 + ms4);
//   //     n5 = (l1 + r5) / (ks1 + ms5);
//   //     n6 = (l1 + r6) / (ks1 + ms6);
//   //     n7 = (l1 + r7) / (ks1 + ms7);
//   //     n8 = (l1 + r8) / (ks1 + ms8);

//   //     d1 = (t1 - l1 - r1) / (k1 * m1) + den_off;
//   //     d2 = (t2 - l1 - r2) / (k1 * m2) + den_off;
//   //     d3 = (t3 - l1 - r3) / (k1 * m3) + den_off;
//   //     d4 = (t4 - l1 - r4) / (k1 * m4) + den_off;
//   //     d5 = (t5 - l1 - r5) / (k1 * m5) + den_off;
//   //     d6 = (t6 - l1 - r6) / (k1 * m6) + den_off;
//   //     d7 = (t7 - l1 - r7) / (k1 * m7) + den_off;
//   //     d8 = (t8 - l1 - r8) / (k1 * m8) + den_off;

//   //     tmp1 = n1 / d1;
//   //     tmp2 = n2 / d2;
//   //     tmp3 = n3 / d3;
//   //     tmp4 = n4 / d4;
//   //     tmp5 = n5 / d5;
//   //     tmp6 = n6 / d6;
//   //     tmp7 = n7 / d7;
//   //     tmp8 = n8 / d8;

//   //     if(tmp1 > maxW){
//   //       maxW = tmp1;
//   //       maxI = ii - lf * gs * gr;
//   //     }
//   //     if(tmp2 > maxW){
//   //       maxW = tmp2;
//   //       maxI = ii - lf * gs * gr + 1;
//   //     }
//   //     if(tmp3 > maxW){
//   //       maxW = tmp3;
//   //       maxI = ii - lf * gs * gr + 2;
//   //     }
//   //     if(tmp4 > maxW){
//   //       maxW = tmp4;
//   //       maxI = ii - lf * gs * gr + 3;
//   //     }
//   //     if(tmp5 > maxW){
//   //       maxW = tmp5;
//   //       maxI = ii - lf * gs * gr + 4;
//   //     }
//   //     if(tmp6 > maxW){
//   //       maxW = tmp6;
//   //       maxI = ii - lf * gs * gr + 5;
//   //     }
//   //     if(tmp7 > maxW){
//   //       maxW = tmp7;
//   //       maxI = ii - lf * gs * gr + 6;
//   //     }
//   //     if(tmp8 > maxW){
//   //       maxW = tmp8;
//   //       maxI = ii - lf * gs * gr + 7;
//   //     }
//   //   }
//   //   for(/*emp*/;j<it_load;j++){
//   //     r1 = rs[ii];       // ii = (wgID * wgSIZE) + (j + (i * it_load) * #wg * wgSIZE) + localID;
//   //     m1 = mss[ii];
//   //     t1 = ts[ii];
//   //     // ii += gs;
//   //     ii += gs * gr;

//   //     ms1 = (m1 * (m1-1)) / 2;
//   //     n1 = (l1 + r1) / (ks1 + ms1);
//   //     d1 = (t1 - l1 - r1) / (k1 * m1) + den_off;
//   //     tmp1 = n1 / d1;

//   //     if(tmp1 > maxW){
//   //       maxW = tmp1;
//   //       maxI = ii - gs * gr;
//   //     }
//   //   }
//   }
//   omega_global[ig] = maxW;
//   index_global[ig] = maxI;
// }

// __kernel void omega16 (
//     __global float *omega_global, __global unsigned int *index_global, __constant float *ls,
//     __constant float *rs, __constant float *ts, __constant int *kss, 
//     __constant int *mss, int gr_load, int it_load, __local float *lls, __local int *lkss
// ) {
//   unsigned int wg = get_group_id(0);
//   unsigned int il = get_local_id(0);
//   unsigned int gs = get_local_size(0);
//   unsigned int ig = get_global_id(0);
//   unsigned int gr = get_num_groups(0);
//   // unsigned int ii = wg * gs * (gr_load * it_load) + il;
//   unsigned int ii = ig;  //wg * gs + il;
//   unsigned int reads = gr_load / gs + 1;
//   unsigned int st = wg * gs * reads;

//   const float den_off = 0.00001f;
//   unsigned int maxI, i, j;

//   float l, r, t, n, d, tmpW, maxW = 0.0f;
//   int k, m, ks, ms;

//   // for(i = 0; i < reads; i++){
//   //   lls[il] = ls[st + il];
//   //   lkss[il] = kss[st + il];
//   //   il += gs;
//   // }

//   lls[il] = ls[st + il];
//   lkss[il] = kss[st + il];

//   // Ensure writes have completed:
//   // barrier(CLK_LOCAL_MEM_FENCE);   // Is stated in optimization guide since work group size is > 64

//   for(i = 0; i < gr_load; i++){
//     l = lls[i];
//     k = lkss[i];
//     ks = (k * (k-1)) / 2;
//     for(j = 0; j < it_load; j++){
//       r = rs[ii];       // ii = (wgID * wgSIZE) + (j + (i * it_load) * #wg * wgSIZE) + localID;
//       m = mss[ii];
//       t = ts[ii];
//       // ii += gs;
//       ii += gs * gr;

//       ms = (m * (m-1)) / 2;
//       n = (l + r) / (ks + ms);
//       d = (t - l - r) / (k * m) + den_off;
//       tmpW = n / d;

//       if(tmpW > maxW){
//         maxW = tmpW;
//         maxI = ii - gs;
//       }
//     }
//   }
//   omega_global[ig] = maxW;
//   index_global[ig] = maxI;

// }

__kernel void omega17 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __constant float *ts, __constant int *km, int outer, int inner
) {
  unsigned int ig = get_global_id(0);
  unsigned int ws = get_local_size(0);
  unsigned int wg = get_group_id(0);

  // unsigned int io = ig & (outer - 1);
  unsigned int io = ig % outer;
  // unsigned int st = wg % (inner / ws) * ws + outer;
  unsigned int st = wg / (outer / ws) * ws + outer;

  // unsigned int st = ((outer/ws) % (inner/ws) == 0 || 
  //                   (inner/ws) % (outer/ws) == 0) ? 
  //                   (ig/inner + wg) % (inner / ws) * ws + outer :
  //                   wg % (inner / ws) * ws + outer;

  const float den_off = 0.00001f;
  unsigned int maxI, i;

  float l, r, t, n, d, tmpW, maxW = 0.0f;
  int k, m, ks, ms;

  l = lr[io];
  k = km[io];
  // ks = (k * (k-1)) / 2;

  #pragma unroll 8
  for(i = st; i < ws + st; i++){
    r = lr[i];
    m = km[i];
    t = ts[i - outer + io * inner];

    // ms = (m * (m-1)) / 2;
    n = (l + r) / ((k * (k-1)) / 2 + (m * (m-1)) / 2);
    d = (t - l - r) / (k * m) + den_off;
    tmpW = n / d;

    if(tmpW > maxW){
      maxW = tmpW;
      maxI = i;
    }
  }
  omega_global[ig] = maxW;
  index_global[ig] = maxI - outer + io * inner;
}

__kernel //__attribute__((reqd_work_group_size(64, 1, 1)))
void omega18 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *lr, 
    __constant float *ts, __constant int *km, int outer, int inner,
    __local float *lrs, __local int *lmss
) {
  unsigned int ig = get_global_id(0);
  unsigned int il = get_local_id(0);
  unsigned int ws = get_local_size(0);
  unsigned int wg = get_group_id(0);

  unsigned int io = ig % outer;
  unsigned int st = (wg / (outer / ws)) * ws + outer;

  const float den_off = 0.00001f;
  unsigned int maxI, i, ii, ip = 0;

  float l, r, t, n, d, tmpW, maxW = 0.0f;
  int k, m, ks, ms;

  l = lr[io];
  k = km[io];
  ks = (k * (k-1)) / 2;

  lrs[il] = lr[st + il];
  lmss[il] = km[st + il];

  // Ensure writes have completed:
  barrier(CLK_LOCAL_MEM_FENCE);   // Is stated in optimization guide since work group size is > 64

  #pragma unroll 8
  for(i = 0; i < ws; i++){
    r = lrs[i];
    m = lmss[i];

    t = ts[i + st - outer + io * inner];

    ms = (m * (m-1)) / 2;
    n = (l + r) / (ks + ms);
    d = (t - l - r) / (k * m) + den_off;
    tmpW = n / d;

    if(tmpW > maxW){
      maxW = tmpW;
      maxI = i;
    }
  }
  omega_global[ig] = maxW;
  index_global[ig] = maxI + st - outer + io * inner;
}

__kernel //__attribute__((reqd_work_group_size(64, 1, 1)))
void omegatest2 (
    __global float *omega_global, __global unsigned int *index_global, __constant float *ls, 
    __constant float *rs, __constant float *ts, __constant int *kss, __constant int *mss,
    int outer, int inner, __local float *lrs, __local int *lmss
) {
  unsigned int ig = get_global_id(0);
  unsigned int gs = get_global_size(0);
  unsigned int il = get_local_id(0);
  unsigned int ws = get_local_size(0);
  unsigned int wg = get_group_id(0);

  unsigned int io = ig % outer;
  unsigned int st = (wg / (inner / ws)) * ws + outer;
  unsigned int stt = ig * ws;
  // unsigned int stt = io * inner + (ig / outer) * ws;
  unsigned int ic = io * inner + st;

  const float den_off = 0.00001f;
  unsigned int maxI, i, ii, ip = 0;

  float l, r, t, n, d, tmpW, maxW = 0.0f;
  int k, m, ks, ms;

  l = ls[io];
  k = kss[io];
  // l = ls[ig];
  // k = kss[ig];
  ks = (k * (k-1)) / 2;

  lrs[il] = ls[st + il];
  lmss[il] = kss[st + il];
  // lrs[il] = rs[ig];
  // lmss[il] = mss[ig];

  // Ensure writes have completed:
  barrier(CLK_LOCAL_MEM_FENCE);   // Is stated in optimization guide since work group size is > 64

  #pragma unroll 8
  for(i = 0; i < ws; i++){
  //   ii = st + ((i + il) & (ws - 1));
  // for(i = st; i < ws + st; i++){
    // ii = (i + il) & (ws - 1);
    r = lrs[i];
    m = lmss[i];
    // r = ls[i];
    // m = kss[i];
    // t = ts[ip + ig];
    // ip += gs;
    t = ts[stt];

    ms = (m * (m-1)) / 2;
    n = (l + r) / (ks + ms);
    d = (t - l - r) / (k * m) + den_off;
    tmpW = n / d;

    if(tmpW > maxW){
      maxW = tmpW;
      // maxI = i;
      maxI = stt;
    }
    stt++;
  }
  omega_global[ig] = maxW;
  index_global[ig] = maxI;
}

__kernel void omega19 (
    __global float *omega, __constant float *rs, __constant float *ts,
    __constant int *m, float l, int k
) {
  const int i = get_global_id(0);

  int ksel2 = k * (k-1) / 2;

  int msel2 = (m[i] * (m[i]-1)) / 2;

  float numerator = (l + rs[i]) / (ksel2 + msel2);

  float denominator = (ts[i] - l - rs[i]) / (k*m[i]) + 0.00001f;

  omega[i] =  numerator / denominator;
}

// On K80 changed __constant to __global is 700us faster?!!?!? Size of __global buffer doesn't matter! Same with HD8750m
__kernel void omega22 (
    __global float *omega, __global unsigned int *index, __global float *LR,
    __global float *TS,__global int *km, int in_cnt, int wi_load
) {
  const unsigned int G_i = get_global_id(0);
  const unsigned int G_s = get_global_size(0);
  unsigned int O_inc = G_s / in_cnt, O_i = G_i / in_cnt + in_cnt, I_i = G_i % in_cnt, G_ic = G_i;

  unsigned int maxI, i;
  float n, d, tmpW, maxW = 0.0f, R_L;
  int k_m, m_k, k_msel2, m_ksel2;

  m_k = km[I_i];
  m_ksel2 = (m_k * (m_k-1)) / 2;

  R_L = LR[I_i];

  #pragma unroll 4
  for(i=0; i<wi_load; i++){
    // O_i = G_ic / in_cnt + in_cnt;
	  // I_i = G_ic % in_cnt;

    k_m = km[O_i];
    k_msel2 = (k_m * (k_m-1)) / 2;

    // m_k = km[I_i];
    // m_ksel2 = (m_k * (m_k-1)) / 2;

    n = (LR[O_i] + R_L) / (k_msel2 + m_ksel2);

    d = (TS[G_ic] - LR[O_i] - R_L) / (k_m * m_k) + 0.00001f;

    tmpW = n / d;

    if(tmpW > maxW){
      maxW = tmpW;
      maxI = G_ic;
    }
    O_i += O_inc;
    G_ic += G_s;
  }

  omega[G_i] = maxW;
  index[G_i] = maxI;
}