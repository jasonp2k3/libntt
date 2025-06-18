#include "ntt/ntt-impl-simd.h"

#define NC 24

static const uint8_t ntt20_fixed_const[NC] = {1, 0, 0, 0, 0, 0,
					      1, 0, 0, 0, 0, 0,
					      1, 0, 0, 0, 0, 0};

extern void X(ntt20_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt20_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, a04, 
	    a05, a06, a07, a08, a09, 
	    a10, a11, a12, a13, a14, 
	    a15, a16, a17, a18, a19;

  {
    sp_simd_t x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 5 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 10 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 15 * istride, idist, vsize);

    a00 = sp_ntt_add_simd(x0, x2, p);
    a01 = sp_ntt_sub_simd(x0, x2, p);
    a02 = sp_ntt_add_simd(x1, x3, p);
    a03 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x0 = sp_simd_gather(in + 4 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 9 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 14 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 19 * istride, idist, vsize);

    a04 = sp_ntt_add_simd(x0, x2, p);
    a05 = sp_ntt_sub_simd(x0, x2, p);
    a06 = sp_ntt_add_simd(x1, x3, p);
    a07 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x3 = sp_simd_gather(in + 3 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 8 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 13 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 18 * istride, idist, vsize);

    a08 = sp_ntt_add_simd(x0, x2, p);
    a09 = sp_ntt_sub_simd(x0, x2, p);
    a10 = sp_ntt_add_simd(x1, x3, p);
    a11 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 7 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 12 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 17 * istride, idist, vsize);

    a12 = sp_ntt_add_simd(x0, x2, p);
    a13 = sp_ntt_sub_simd(x0, x2, p);
    a14 = sp_ntt_add_simd(x1, x3, p);
    a15 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 6 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 11 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 16 * istride, idist, vsize);

    a16 = sp_ntt_add_simd(x0, x2, p);
    a17 = sp_ntt_sub_simd(x0, x2, p);
    a18 = sp_ntt_add_simd(x1, x3, p);
    a19 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a04;
    b4 = a08;
    b2 = a12;
    b3 = a16;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a00 = b0;
    a04 = b4;
    a08 = b3;
    a12 = b1;
    a16 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a05;
    b4 = a09;
    b2 = a13;
    b3 = a17;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a01 = b0;
    a05 = b4;
    a09 = b3;
    a13 = b1;
    a17 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a06;
    b4 = a10;
    b2 = a14;
    b3 = a18;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a02 = b0;
    a06 = b4;
    a10 = b3;
    a14 = b1;
    a18 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a07;
    b4 = a11;
    b2 = a15;
    b3 = a19;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b0 = sp_ntt_mul_simd(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a03 = b0;
    a07 = b4;
    a11 = b3;
    a15 = b1;
    a19 = b2;
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a00, a02, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a00, a02, p);
    sp_simd_t p2 = sp_ntt_add_simd(a01, a03, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a01, a03, p);

    sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 5 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 10 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 15 * ostride, odist, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a04, a06, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a04, a06, p);
    sp_simd_t p2 = sp_ntt_add_simd(a05, a07, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a05, a07, p);

    sp_simd_scatter(p2, out + 1 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 6 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 11 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 16 * ostride, odist, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a08, a10, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a08, a10, p);
    sp_simd_t p2 = sp_ntt_add_simd(a09, a11, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a09, a11, p);

    sp_simd_scatter(p1, out + 2 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 7 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 12 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 17 * ostride, odist, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a12, a14, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a12, a14, p);
    sp_simd_t p2 = sp_ntt_add_simd(a13, a15, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a13, a15, p);

    sp_simd_scatter(p3, out + 3 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 8 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 13 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 18 * ostride, odist, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a16, a18, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a16, a18, p);
    sp_simd_t p2 = sp_ntt_add_simd(a17, a19, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a17, a19, p);

    sp_simd_scatter(p0, out + 4 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 9 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 14 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 19 * ostride, odist, vsize);
  }
}

static void
ntt20_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, a04, 
	    a05, a06, a07, a08, a09, 
	    a10, a11, a12, a13, a14, 
	    a15, a16, a17, a18, a19;

  {
    sp_simd_t x0 = sp_simd_load(in + 0 * istride);
    sp_simd_t x1 = sp_simd_load(in + 5 * istride);
    sp_simd_t x2 = sp_simd_load(in + 10 * istride);
    sp_simd_t x3 = sp_simd_load(in + 15 * istride);

    a00 = sp_ntt_add_simd0(x0, x2, p);
    a01 = sp_ntt_sub_simd0(x0, x2, p);
    a02 = sp_ntt_add_simd0(x1, x3, p);
    a03 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x0 = sp_simd_load(in + 4 * istride);
    sp_simd_t x1 = sp_simd_load(in + 9 * istride);
    sp_simd_t x2 = sp_simd_load(in + 14 * istride);
    sp_simd_t x3 = sp_simd_load(in + 19 * istride);

    a04 = sp_ntt_add_simd0(x0, x2, p);
    a05 = sp_ntt_sub_simd0(x0, x2, p);
    a06 = sp_ntt_add_simd0(x1, x3, p);
    a07 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x3 = sp_simd_load(in + 3 * istride);
    sp_simd_t x0 = sp_simd_load(in + 8 * istride);
    sp_simd_t x1 = sp_simd_load(in + 13 * istride);
    sp_simd_t x2 = sp_simd_load(in + 18 * istride);

    a08 = sp_ntt_add_simd0(x0, x2, p);
    a09 = sp_ntt_sub_simd0(x0, x2, p);
    a10 = sp_ntt_add_simd0(x1, x3, p);
    a11 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x2 = sp_simd_load(in + 2 * istride);
    sp_simd_t x3 = sp_simd_load(in + 7 * istride);
    sp_simd_t x0 = sp_simd_load(in + 12 * istride);
    sp_simd_t x1 = sp_simd_load(in + 17 * istride);

    a12 = sp_ntt_add_simd0(x0, x2, p);
    a13 = sp_ntt_sub_simd0(x0, x2, p);
    a14 = sp_ntt_add_simd0(x1, x3, p);
    a15 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x1 = sp_simd_load(in + 1 * istride);
    sp_simd_t x2 = sp_simd_load(in + 6 * istride);
    sp_simd_t x3 = sp_simd_load(in + 11 * istride);
    sp_simd_t x0 = sp_simd_load(in + 16 * istride);

    a16 = sp_ntt_add_simd0(x0, x2, p);
    a17 = sp_ntt_sub_simd0(x0, x2, p);
    a18 = sp_ntt_add_simd0(x1, x3, p);
    a19 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a04;
    b4 = a08;
    b2 = a12;
    b3 = a16;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c + 2*1, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*2, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*3, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*4, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*5, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a00 = b0;
    a04 = b4;
    a08 = b3;
    a12 = b1;
    a16 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a05;
    b4 = a09;
    b2 = a13;
    b3 = a17;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c + 2*7, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*8, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*9, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*10, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*11, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a01 = b0;
    a05 = b4;
    a09 = b3;
    a13 = b1;
    a17 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a06;
    b4 = a10;
    b2 = a14;
    b3 = a18;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c + 2*13, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*14, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*15, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*16, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*17, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a02 = b0;
    a06 = b4;
    a10 = b3;
    a14 = b1;
    a18 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a07;
    b4 = a11;
    b2 = a15;
    b3 = a19;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b0 = sp_ntt_mul_simd0(b0, c + 2*18, p);
    b1 = sp_ntt_mul_simd0(b1, c + 2*19, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*20, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*21, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*22, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*23, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a03 = b0;
    a07 = b4;
    a11 = b3;
    a15 = b1;
    a19 = b2;
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a00, a02, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a00, a02, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a01, a03, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a01, a03, p);

    sp_simd_store(p0, out + 0 * ostride);
    sp_simd_store(p2, out + 5 * ostride);
    sp_simd_store(p1, out + 10 * ostride);
    sp_simd_store(p3, out + 15 * ostride);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a04, a06, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a04, a06, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a05, a07, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a05, a07, p);

    sp_simd_store(p2, out + 1 * ostride);
    sp_simd_store(p1, out + 6 * ostride);
    sp_simd_store(p3, out + 11 * ostride);
    sp_simd_store(p0, out + 16 * ostride);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a08, a10, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a08, a10, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a09, a11, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a09, a11, p);

    sp_simd_store(p1, out + 2 * ostride);
    sp_simd_store(p3, out + 7 * ostride);
    sp_simd_store(p0, out + 12 * ostride);
    sp_simd_store(p2, out + 17 * ostride);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a12, a14, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a12, a14, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a13, a15, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a13, a15, p);

    sp_simd_store(p3, out + 3 * ostride);
    sp_simd_store(p0, out + 8 * ostride);
    sp_simd_store(p2, out + 13 * ostride);
    sp_simd_store(p1, out + 18 * ostride);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a16, a18, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a16, a18, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a17, a19, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a17, a19, p);

    sp_simd_store(p0, out + 4 * ostride);
    sp_simd_store(p2, out + 9 * ostride);
    sp_simd_store(p1, out + 14 * ostride);
    sp_simd_store(p3, out + 19 * ostride);
  }
}


static void
ntt20_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, a04, 
	    a05, a06, a07, a08, a09, 
	    a10, a11, a12, a13, a14, 
	    a15, a16, a17, a18, a19;

  {
    sp_simd_t x0 = sp_simd_gather(in + 0 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 5 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 10 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 15 * istride, idist, vsize);

    a00 = sp_ntt_add_simd(x0, x2, p);
    a01 = sp_ntt_sub_simd(x0, x2, p);
    a02 = sp_ntt_add_simd(x1, x3, p);
    a03 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x0 = sp_simd_gather(in + 4 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 9 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 14 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 19 * istride, idist, vsize);

    a04 = sp_ntt_add_simd(x0, x2, p);
    a05 = sp_ntt_sub_simd(x0, x2, p);
    a06 = sp_ntt_add_simd(x1, x3, p);
    a07 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x3 = sp_simd_gather(in + 3 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 8 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 13 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 18 * istride, idist, vsize);

    a08 = sp_ntt_add_simd(x0, x2, p);
    a09 = sp_ntt_sub_simd(x0, x2, p);
    a10 = sp_ntt_add_simd(x1, x3, p);
    a11 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x2 = sp_simd_gather(in + 2 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 7 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 12 * istride, idist, vsize);
    sp_simd_t x1 = sp_simd_gather(in + 17 * istride, idist, vsize);

    a12 = sp_ntt_add_simd(x0, x2, p);
    a13 = sp_ntt_sub_simd(x0, x2, p);
    a14 = sp_ntt_add_simd(x1, x3, p);
    a15 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x1 = sp_simd_gather(in + 1 * istride, idist, vsize);
    sp_simd_t x2 = sp_simd_gather(in + 6 * istride, idist, vsize);
    sp_simd_t x3 = sp_simd_gather(in + 11 * istride, idist, vsize);
    sp_simd_t x0 = sp_simd_gather(in + 16 * istride, idist, vsize);

    a16 = sp_ntt_add_simd(x0, x2, p);
    a17 = sp_ntt_sub_simd(x0, x2, p);
    a18 = sp_ntt_add_simd(x1, x3, p);
    a19 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a04;
    b4 = a08;
    b2 = a12;
    b3 = a16;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a00 = b0;
    a04 = b4;
    a08 = b3;
    a12 = b1;
    a16 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a05;
    b4 = a09;
    b2 = a13;
    b3 = a17;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a01 = b0;
    a05 = b4;
    a09 = b3;
    a13 = b1;
    a17 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a06;
    b4 = a10;
    b2 = a14;
    b3 = a18;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a02 = b0;
    a06 = b4;
    a10 = b3;
    a14 = b1;
    a18 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a07;
    b4 = a11;
    b2 = a15;
    b3 = a19;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b0 = sp_ntt_mul_simd(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a03 = b0;
    a07 = b4;
    a11 = b3;
    a15 = b1;
    a19 = b2;
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a00, a02, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd(a00, a02, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd(a01, a03, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd(a01, a03, p);

    p2 = sp_ntt_twiddle_mul_simd(p2, w + 8, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 18, p);
    p3 = sp_ntt_twiddle_mul_simd(p3, w + 28, p);

    sp_simd_scatter(p0, out + 0 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 5 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 10 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 15 * ostride, odist, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_partial_simd(a04, a06, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd(a04, a06, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd(a05, a07, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd(a05, a07, p);

    p2 = sp_ntt_twiddle_mul_simd(p2, w + 0, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 10, p);
    p3 = sp_ntt_twiddle_mul_simd(p3, w + 20, p);
    p0 = sp_ntt_twiddle_mul_simd(p0, w + 30, p);

    sp_simd_scatter(p2, out + 1 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 6 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 11 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 16 * ostride, odist, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_partial_simd(a08, a10, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd(a08, a10, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd(a09, a11, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd(a09, a11, p);

    p1 = sp_ntt_twiddle_mul_simd(p1, w + 2, p);
    p3 = sp_ntt_twiddle_mul_simd(p3, w + 12, p);
    p0 = sp_ntt_twiddle_mul_simd(p0, w + 22, p);
    p2 = sp_ntt_twiddle_mul_simd(p2, w + 32, p);

    sp_simd_scatter(p1, out + 2 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 7 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 12 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 17 * ostride, odist, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_partial_simd(a12, a14, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd(a12, a14, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd(a13, a15, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd(a13, a15, p);

    p3 = sp_ntt_twiddle_mul_simd(p3, w + 4, p);
    p0 = sp_ntt_twiddle_mul_simd(p0, w + 14, p);
    p2 = sp_ntt_twiddle_mul_simd(p2, w + 24, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 34, p);

    sp_simd_scatter(p3, out + 3 * ostride, odist, vsize);
    sp_simd_scatter(p0, out + 8 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 13 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 18 * ostride, odist, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_partial_simd(a16, a18, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd(a16, a18, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd(a17, a19, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd(a17, a19, p);

    p0 = sp_ntt_twiddle_mul_simd(p0, w + 6, p);
    p2 = sp_ntt_twiddle_mul_simd(p2, w + 16, p);
    p1 = sp_ntt_twiddle_mul_simd(p1, w + 26, p);
    p3 = sp_ntt_twiddle_mul_simd(p3, w + 36, p);

    sp_simd_scatter(p0, out + 4 * ostride, odist, vsize);
    sp_simd_scatter(p2, out + 9 * ostride, odist, vsize);
    sp_simd_scatter(p1, out + 14 * ostride, odist, vsize);
    sp_simd_scatter(p3, out + 19 * ostride, odist, vsize);
  }
}

static void
ntt20_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, a04, 
	    a05, a06, a07, a08, a09, 
	    a10, a11, a12, a13, a14, 
	    a15, a16, a17, a18, a19;

  {
    sp_simd_t x0 = sp_simd_load(in + 0 * istride);
    sp_simd_t x1 = sp_simd_load(in + 5 * istride);
    sp_simd_t x2 = sp_simd_load(in + 10 * istride);
    sp_simd_t x3 = sp_simd_load(in + 15 * istride);

    a00 = sp_ntt_add_simd0(x0, x2, p);
    a01 = sp_ntt_sub_simd0(x0, x2, p);
    a02 = sp_ntt_add_simd0(x1, x3, p);
    a03 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x0 = sp_simd_load(in + 4 * istride);
    sp_simd_t x1 = sp_simd_load(in + 9 * istride);
    sp_simd_t x2 = sp_simd_load(in + 14 * istride);
    sp_simd_t x3 = sp_simd_load(in + 19 * istride);

    a04 = sp_ntt_add_simd0(x0, x2, p);
    a05 = sp_ntt_sub_simd0(x0, x2, p);
    a06 = sp_ntt_add_simd0(x1, x3, p);
    a07 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x3 = sp_simd_load(in + 3 * istride);
    sp_simd_t x0 = sp_simd_load(in + 8 * istride);
    sp_simd_t x1 = sp_simd_load(in + 13 * istride);
    sp_simd_t x2 = sp_simd_load(in + 18 * istride);

    a08 = sp_ntt_add_simd0(x0, x2, p);
    a09 = sp_ntt_sub_simd0(x0, x2, p);
    a10 = sp_ntt_add_simd0(x1, x3, p);
    a11 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x2 = sp_simd_load(in + 2 * istride);
    sp_simd_t x3 = sp_simd_load(in + 7 * istride);
    sp_simd_t x0 = sp_simd_load(in + 12 * istride);
    sp_simd_t x1 = sp_simd_load(in + 17 * istride);

    a12 = sp_ntt_add_simd0(x0, x2, p);
    a13 = sp_ntt_sub_simd0(x0, x2, p);
    a14 = sp_ntt_add_simd0(x1, x3, p);
    a15 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x1 = sp_simd_load(in + 1 * istride);
    sp_simd_t x2 = sp_simd_load(in + 6 * istride);
    sp_simd_t x3 = sp_simd_load(in + 11 * istride);
    sp_simd_t x0 = sp_simd_load(in + 16 * istride);

    a16 = sp_ntt_add_simd0(x0, x2, p);
    a17 = sp_ntt_sub_simd0(x0, x2, p);
    a18 = sp_ntt_add_simd0(x1, x3, p);
    a19 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a04;
    b4 = a08;
    b2 = a12;
    b3 = a16;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c + 2*1, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*2, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*3, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*4, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*5, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a00 = b0;
    a04 = b4;
    a08 = b3;
    a12 = b1;
    a16 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a05;
    b4 = a09;
    b2 = a13;
    b3 = a17;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c + 2*7, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*8, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*9, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*10, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*11, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a01 = b0;
    a05 = b4;
    a09 = b3;
    a13 = b1;
    a17 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a06;
    b4 = a10;
    b2 = a14;
    b3 = a18;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c + 2*13, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*14, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*15, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*16, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*17, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a02 = b0;
    a06 = b4;
    a10 = b3;
    a14 = b1;
    a18 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a07;
    b4 = a11;
    b2 = a15;
    b3 = a19;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b0 = sp_ntt_mul_simd0(b0, c + 2*18, p);
    b1 = sp_ntt_mul_simd0(b1, c + 2*19, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*20, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*21, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*22, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*23, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a03 = b0;
    a07 = b4;
    a11 = b3;
    a15 = b1;
    a19 = b2;
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a00, a02, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd0(a00, a02, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd0(a01, a03, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd0(a01, a03, p);

    p2 = sp_ntt_twiddle_mul_simd0(p2, w + 8, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w + 18, p);
    p3 = sp_ntt_twiddle_mul_simd0(p3, w + 28, p);

    sp_simd_store(p0, out + 0 * ostride);
    sp_simd_store(p2, out + 5 * ostride);
    sp_simd_store(p1, out + 10 * ostride);
    sp_simd_store(p3, out + 15 * ostride);
  }
  {
    sp_simd_t p0 = sp_ntt_add_partial_simd0(a04, a06, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd0(a04, a06, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd0(a05, a07, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd0(a05, a07, p);

    p2 = sp_ntt_twiddle_mul_simd0(p2, w + 0, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w + 10, p);
    p3 = sp_ntt_twiddle_mul_simd0(p3, w + 20, p);
    p0 = sp_ntt_twiddle_mul_simd0(p0, w + 30, p);

    sp_simd_store(p2, out + 1 * ostride);
    sp_simd_store(p1, out + 6 * ostride);
    sp_simd_store(p3, out + 11 * ostride);
    sp_simd_store(p0, out + 16 * ostride);
  }
  {
    sp_simd_t p0 = sp_ntt_add_partial_simd0(a08, a10, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd0(a08, a10, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd0(a09, a11, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd0(a09, a11, p);

    p1 = sp_ntt_twiddle_mul_simd0(p1, w + 2, p);
    p3 = sp_ntt_twiddle_mul_simd0(p3, w + 12, p);
    p0 = sp_ntt_twiddle_mul_simd0(p0, w + 22, p);
    p2 = sp_ntt_twiddle_mul_simd0(p2, w + 32, p);

    sp_simd_store(p1, out + 2 * ostride);
    sp_simd_store(p3, out + 7 * ostride);
    sp_simd_store(p0, out + 12 * ostride);
    sp_simd_store(p2, out + 17 * ostride);
  }
  {
    sp_simd_t p0 = sp_ntt_add_partial_simd0(a12, a14, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd0(a12, a14, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd0(a13, a15, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd0(a13, a15, p);

    p3 = sp_ntt_twiddle_mul_simd0(p3, w + 4, p);
    p0 = sp_ntt_twiddle_mul_simd0(p0, w + 14, p);
    p2 = sp_ntt_twiddle_mul_simd0(p2, w + 24, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w + 34, p);

    sp_simd_store(p3, out + 3 * ostride);
    sp_simd_store(p0, out + 8 * ostride);
    sp_simd_store(p2, out + 13 * ostride);
    sp_simd_store(p1, out + 18 * ostride);
  }
  {
    sp_simd_t p0 = sp_ntt_add_partial_simd0(a16, a18, p);
    sp_simd_t p1 = sp_ntt_sub_partial_simd0(a16, a18, p);
    sp_simd_t p2 = sp_ntt_add_partial_simd0(a17, a19, p);
    sp_simd_t p3 = sp_ntt_sub_partial_simd0(a17, a19, p);

    p0 = sp_ntt_twiddle_mul_simd0(p0, w + 6, p);
    p2 = sp_ntt_twiddle_mul_simd0(p2, w + 16, p);
    p1 = sp_ntt_twiddle_mul_simd0(p1, w + 26, p);
    p3 = sp_ntt_twiddle_mul_simd0(p3, w + 36, p);

    sp_simd_store(p0, out + 4 * ostride);
    sp_simd_store(p2, out + 9 * ostride);
    sp_simd_store(p1, out + 14 * ostride);
    sp_simd_store(p3, out + 19 * ostride);
  }
}


static void
ntt20_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j00, j01, j02, j03, j04,
             j05, j06, j07, j08, j09,
	     j10, j11, j12, j13, j14,
	     j15, j16, j17, j18, j19;

  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14,
            a15, a16, a17, a18, a19;

  j00 = start;
  j01 = sp_array_inc(j00, inc, n);
  j02 = sp_array_inc(j00, 2 * inc, n);
  j03 = sp_array_inc(j00, 3 * inc, n);
  j04 = sp_array_inc(j00, 4 * inc, n);
  j05 = sp_array_inc(j00, 5 * inc, n);
  j06 = sp_array_inc(j00, 6 * inc, n);
  j07 = sp_array_inc(j00, 7 * inc, n);
  j08 = sp_array_inc(j00, 8 * inc, n);
  j09 = sp_array_inc(j00, 9 * inc, n);
  j10 = sp_array_inc(j00, 10 * inc, n);
  j11 = sp_array_inc(j00, 11 * inc, n);
  j12 = sp_array_inc(j00, 12 * inc, n);
  j13 = sp_array_inc(j00, 13 * inc, n);
  j14 = sp_array_inc(j00, 14 * inc, n);
  j15 = sp_array_inc(j00, 15 * inc, n);
  j16 = sp_array_inc(j00, 16 * inc, n);
  j17 = sp_array_inc(j00, 17 * inc, n);
  j18 = sp_array_inc(j00, 18 * inc, n);
  j19 = sp_array_inc(j00, 19 * inc, n);

  {
    sp_simd_t x0 = sp_simd_pfa_gather(x, j00, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j05, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j10, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j15, inc2, n, vsize);

    a00 = sp_ntt_add_simd(x0, x2, p);
    a01 = sp_ntt_sub_simd(x0, x2, p);
    a02 = sp_ntt_add_simd(x1, x3, p);
    a03 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x0 = sp_simd_pfa_gather(x, j04, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j09, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j14, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j19, inc2, n, vsize);

    a04 = sp_ntt_add_simd(x0, x2, p);
    a05 = sp_ntt_sub_simd(x0, x2, p);
    a06 = sp_ntt_add_simd(x1, x3, p);
    a07 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x3 = sp_simd_pfa_gather(x, j03, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j08, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j13, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j18, inc2, n, vsize);

    a08 = sp_ntt_add_simd(x0, x2, p);
    a09 = sp_ntt_sub_simd(x0, x2, p);
    a10 = sp_ntt_add_simd(x1, x3, p);
    a11 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x2 = sp_simd_pfa_gather(x, j02, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j07, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j12, inc2, n, vsize);
    sp_simd_t x1 = sp_simd_pfa_gather(x, j17, inc2, n, vsize);

    a12 = sp_ntt_add_simd(x0, x2, p);
    a13 = sp_ntt_sub_simd(x0, x2, p);
    a14 = sp_ntt_add_simd(x1, x3, p);
    a15 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t x1 = sp_simd_pfa_gather(x, j01, inc2, n, vsize);
    sp_simd_t x2 = sp_simd_pfa_gather(x, j06, inc2, n, vsize);
    sp_simd_t x3 = sp_simd_pfa_gather(x, j11, inc2, n, vsize);
    sp_simd_t x0 = sp_simd_pfa_gather(x, j16, inc2, n, vsize);

    a16 = sp_ntt_add_simd(x0, x2, p);
    a17 = sp_ntt_sub_simd(x0, x2, p);
    a18 = sp_ntt_add_simd(x1, x3, p);
    a19 = sp_ntt_sub_simd(x1, x3, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a04;
    b4 = a08;
    b2 = a12;
    b3 = a16;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a00 = b0;
    a04 = b4;
    a08 = b3;
    a12 = b1;
    a16 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a05;
    b4 = a09;
    b2 = a13;
    b3 = a17;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a01 = b0;
    a05 = b4;
    a09 = b3;
    a13 = b1;
    a17 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a06;
    b4 = a10;
    b2 = a14;
    b3 = a18;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b1 = sp_ntt_mul_simd(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a02 = b0;
    a06 = b4;
    a10 = b3;
    a14 = b1;
    a18 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a07;
    b4 = a11;
    b2 = a15;
    b3 = a19;

    c1 = sp_ntt_add_simd(b1, b3, p);
    c3 = sp_ntt_sub_simd(b1, b3, p);
    c2 = sp_ntt_add_simd(b2, b4, p);
    c4 = sp_ntt_sub_simd(b2, b4, p);

    b1 = sp_ntt_add_simd(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd(c3, c4, p);

    b0 = sp_ntt_add_simd(b0, b1, p);

    b0 = sp_ntt_mul_simd(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul_simd(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul_simd(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul_simd(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul_simd(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul_simd(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add_simd(b0, b1, p);

    c1 = sp_ntt_add_simd(b1, b2, p);
    c2 = sp_ntt_sub_simd(b1, b2, p);
    c3 = sp_ntt_add_simd(b3, b5, p);
    c4 = sp_ntt_add_simd(b4, b5, p);

    b1 = sp_ntt_add_simd(c1, c3, p);
    b2 = sp_ntt_add_simd(c2, c4, p);
    b3 = sp_ntt_sub_simd(c1, c3, p);
    b4 = sp_ntt_sub_simd(c2, c4, p);

    a03 = b0;
    a07 = b4;
    a11 = b3;
    a15 = b1;
    a19 = b2;
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a00, a02, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a00, a02, p);
    sp_simd_t p2 = sp_ntt_add_simd(a01, a03, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a01, a03, p);

    sp_simd_pfa_scatter(p0, x, j00, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j05, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j10, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j15, inc2, n, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a04, a06, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a04, a06, p);
    sp_simd_t p2 = sp_ntt_add_simd(a05, a07, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a05, a07, p);

    sp_simd_pfa_scatter(p2, x, j01, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j06, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j11, inc2, n, vsize);
    sp_simd_pfa_scatter(p0, x, j16, inc2, n, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a08, a10, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a08, a10, p);
    sp_simd_t p2 = sp_ntt_add_simd(a09, a11, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a09, a11, p);

    sp_simd_pfa_scatter(p1, x, j02, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j07, inc2, n, vsize);
    sp_simd_pfa_scatter(p0, x, j12, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j17, inc2, n, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a12, a14, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a12, a14, p);
    sp_simd_t p2 = sp_ntt_add_simd(a13, a15, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a13, a15, p);

    sp_simd_pfa_scatter(p3, x, j03, inc2, n, vsize);
    sp_simd_pfa_scatter(p0, x, j08, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j13, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j18, inc2, n, vsize);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd(a16, a18, p);
    sp_simd_t p1 = sp_ntt_sub_simd(a16, a18, p);
    sp_simd_t p2 = sp_ntt_add_simd(a17, a19, p);
    sp_simd_t p3 = sp_ntt_sub_simd(a17, a19, p);

    sp_simd_pfa_scatter(p0, x, j04, inc2, n, vsize);
    sp_simd_pfa_scatter(p2, x, j09, inc2, n, vsize);
    sp_simd_pfa_scatter(p1, x, j14, inc2, n, vsize);
    sp_simd_pfa_scatter(p3, x, j19, inc2, n, vsize);
  }
}

static void
ntt20_pfa_run_core_simd_interleaved(
	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j00, j01, j02, j03, j04,
             j05, j06, j07, j08, j09,
	     j10, j11, j12, j13, j14,
	     j15, j16, j17, j18, j19;

  sp_simd_t a00, a01, a02, a03, a04, 
            a05, a06, a07, a08, a09, 
            a10, a11, a12, a13, a14,
            a15, a16, a17, a18, a19;

  j00 = start;
  j01 = sp_array_inc(j00, inc, n);
  j02 = sp_array_inc(j00, 2 * inc, n);
  j03 = sp_array_inc(j00, 3 * inc, n);
  j04 = sp_array_inc(j00, 4 * inc, n);
  j05 = sp_array_inc(j00, 5 * inc, n);
  j06 = sp_array_inc(j00, 6 * inc, n);
  j07 = sp_array_inc(j00, 7 * inc, n);
  j08 = sp_array_inc(j00, 8 * inc, n);
  j09 = sp_array_inc(j00, 9 * inc, n);
  j10 = sp_array_inc(j00, 10 * inc, n);
  j11 = sp_array_inc(j00, 11 * inc, n);
  j12 = sp_array_inc(j00, 12 * inc, n);
  j13 = sp_array_inc(j00, 13 * inc, n);
  j14 = sp_array_inc(j00, 14 * inc, n);
  j15 = sp_array_inc(j00, 15 * inc, n);
  j16 = sp_array_inc(j00, 16 * inc, n);
  j17 = sp_array_inc(j00, 17 * inc, n);
  j18 = sp_array_inc(j00, 18 * inc, n);
  j19 = sp_array_inc(j00, 19 * inc, n);

  {
    sp_simd_t x0 = sp_simd_load(x + j00);
    sp_simd_t x1 = sp_simd_load(x + j05);
    sp_simd_t x2 = sp_simd_load(x + j10);
    sp_simd_t x3 = sp_simd_load(x + j15);

    a00 = sp_ntt_add_simd0(x0, x2, p);
    a01 = sp_ntt_sub_simd0(x0, x2, p);
    a02 = sp_ntt_add_simd0(x1, x3, p);
    a03 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x0 = sp_simd_load(x + j04);
    sp_simd_t x1 = sp_simd_load(x + j09);
    sp_simd_t x2 = sp_simd_load(x + j14);
    sp_simd_t x3 = sp_simd_load(x + j19);

    a04 = sp_ntt_add_simd0(x0, x2, p);
    a05 = sp_ntt_sub_simd0(x0, x2, p);
    a06 = sp_ntt_add_simd0(x1, x3, p);
    a07 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x3 = sp_simd_load(x + j03);
    sp_simd_t x0 = sp_simd_load(x + j08);
    sp_simd_t x1 = sp_simd_load(x + j13);
    sp_simd_t x2 = sp_simd_load(x + j18);

    a08 = sp_ntt_add_simd0(x0, x2, p);
    a09 = sp_ntt_sub_simd0(x0, x2, p);
    a10 = sp_ntt_add_simd0(x1, x3, p);
    a11 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x2 = sp_simd_load(x + j02);
    sp_simd_t x3 = sp_simd_load(x + j07);
    sp_simd_t x0 = sp_simd_load(x + j12);
    sp_simd_t x1 = sp_simd_load(x + j17);

    a12 = sp_ntt_add_simd0(x0, x2, p);
    a13 = sp_ntt_sub_simd0(x0, x2, p);
    a14 = sp_ntt_add_simd0(x1, x3, p);
    a15 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t x1 = sp_simd_load(x + j01);
    sp_simd_t x2 = sp_simd_load(x + j06);
    sp_simd_t x3 = sp_simd_load(x + j11);
    sp_simd_t x0 = sp_simd_load(x + j16);

    a16 = sp_ntt_add_simd0(x0, x2, p);
    a17 = sp_ntt_sub_simd0(x0, x2, p);
    a18 = sp_ntt_add_simd0(x1, x3, p);
    a19 = sp_ntt_sub_simd0(x1, x3, p);
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a04;
    b4 = a08;
    b2 = a12;
    b3 = a16;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c + 2*1, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*2, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*3, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*4, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*5, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a00 = b0;
    a04 = b4;
    a08 = b3;
    a12 = b1;
    a16 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a05;
    b4 = a09;
    b2 = a13;
    b3 = a17;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c + 2*7, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*8, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*9, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*10, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*11, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a01 = b0;
    a05 = b4;
    a09 = b3;
    a13 = b1;
    a17 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a06;
    b4 = a10;
    b2 = a14;
    b3 = a18;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b1 = sp_ntt_mul_simd0(b1, c + 2*13, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*14, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*15, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*16, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*17, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a02 = b0;
    a06 = b4;
    a10 = b3;
    a14 = b1;
    a18 = b2;
  }
  {
    sp_simd_t b0, b1, b2, b3, b4, b5;
    sp_simd_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a07;
    b4 = a11;
    b2 = a15;
    b3 = a19;

    c1 = sp_ntt_add_simd0(b1, b3, p);
    c3 = sp_ntt_sub_simd0(b1, b3, p);
    c2 = sp_ntt_add_simd0(b2, b4, p);
    c4 = sp_ntt_sub_simd0(b2, b4, p);

    b1 = sp_ntt_add_simd0(c1, c2, p);
    b2 = sp_ntt_sub_partial_simd0(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial_simd0(c3, c4, p);

    b0 = sp_ntt_add_simd0(b0, b1, p);

    b0 = sp_ntt_mul_simd0(b0, c + 2*18, p);
    b1 = sp_ntt_mul_simd0(b1, c + 2*19, p);
    b2 = sp_ntt_mul_simd0(b2, c + 2*20, p);
    b3 = sp_ntt_mul_simd0(b3, c + 2*21, p);
    b4 = sp_ntt_mul_simd0(b4, c + 2*22, p);
    b5 = sp_ntt_mul_simd0(b5, c + 2*23, p);

    b1 = sp_ntt_add_simd0(b0, b1, p);

    c1 = sp_ntt_add_simd0(b1, b2, p);
    c2 = sp_ntt_sub_simd0(b1, b2, p);
    c3 = sp_ntt_add_simd0(b3, b5, p);
    c4 = sp_ntt_add_simd0(b4, b5, p);

    b1 = sp_ntt_add_simd0(c1, c3, p);
    b2 = sp_ntt_add_simd0(c2, c4, p);
    b3 = sp_ntt_sub_simd0(c1, c3, p);
    b4 = sp_ntt_sub_simd0(c2, c4, p);

    a03 = b0;
    a07 = b4;
    a11 = b3;
    a15 = b1;
    a19 = b2;
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a00, a02, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a00, a02, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a01, a03, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a01, a03, p);

    sp_simd_store(p0, x + j00);
    sp_simd_store(p2, x + j05);
    sp_simd_store(p1, x + j10);
    sp_simd_store(p3, x + j15);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a04, a06, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a04, a06, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a05, a07, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a05, a07, p);

    sp_simd_store(p2, x + j01);
    sp_simd_store(p1, x + j06);
    sp_simd_store(p3, x + j11);
    sp_simd_store(p0, x + j16);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a08, a10, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a08, a10, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a09, a11, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a09, a11, p);

    sp_simd_store(p1, x + j02);
    sp_simd_store(p3, x + j07);
    sp_simd_store(p0, x + j12);
    sp_simd_store(p2, x + j17);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a12, a14, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a12, a14, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a13, a15, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a13, a15, p);

    sp_simd_store(p3, x + j03);
    sp_simd_store(p0, x + j08);
    sp_simd_store(p2, x + j13);
    sp_simd_store(p1, x + j18);
  }
  {
    sp_simd_t p0 = sp_ntt_add_simd0(a16, a18, p);
    sp_simd_t p1 = sp_ntt_sub_simd0(a16, a18, p);
    sp_simd_t p2 = sp_ntt_add_simd0(a17, a19, p);
    sp_simd_t p3 = sp_ntt_sub_simd0(a17, a19, p);

    sp_simd_store(p0, x + j04);
    sp_simd_store(p2, x + j09);
    sp_simd_store(p1, x + j14);
    sp_simd_store(p3, x + j19);
  }
}

DECLARE_CORE_ROUTINES_SIMD(20)
