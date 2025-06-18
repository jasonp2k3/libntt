#include "ntt/ntt-impl-scalar.h"

#define NC 24

static const uint8_t ntt20_fixed_const[NC] = {1, 0, 0, 0, 0, 0,
					      1, 0, 0, 0, 0, 0,
					      1, 0, 0, 0, 0, 0};

void
X(ntt20_init)(spv_t out, sp_t p, sp_t d,
	  sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt20_config), out, p, d, primroot, order, perm);
}

static void 
ntt20_run_core(spv_t in, spv_size_t istride,
                spv_t out, spv_size_t ostride,
                sp_t p, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, a04, 
       a05, a06, a07, a08, a09, 
       a10, a11, a12, a13, a14, 
       a15, a16, a17, a18, a19;

  {
    sp_t x0 = in[0 * istride];
    sp_t x1 = in[5 * istride];
    sp_t x2 = in[10 * istride];
    sp_t x3 = in[15 * istride];

    a00 = sp_ntt_add(x0, x2, p);
    a01 = sp_ntt_sub(x0, x2, p);
    a02 = sp_ntt_add(x1, x3, p);
    a03 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x0 = in[4 * istride];
    sp_t x1 = in[9 * istride];
    sp_t x2 = in[14 * istride];
    sp_t x3 = in[19 * istride];

    a04 = sp_ntt_add(x0, x2, p);
    a05 = sp_ntt_sub(x0, x2, p);
    a06 = sp_ntt_add(x1, x3, p);
    a07 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x3 = in[3 * istride];
    sp_t x0 = in[8 * istride];
    sp_t x1 = in[13 * istride];
    sp_t x2 = in[18 * istride];

    a08 = sp_ntt_add(x0, x2, p);
    a09 = sp_ntt_sub(x0, x2, p);
    a10 = sp_ntt_add(x1, x3, p);
    a11 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x2 = in[2 * istride];
    sp_t x3 = in[7 * istride];
    sp_t x0 = in[12 * istride];
    sp_t x1 = in[17 * istride];

    a12 = sp_ntt_add(x0, x2, p);
    a13 = sp_ntt_sub(x0, x2, p);
    a14 = sp_ntt_add(x1, x3, p);
    a15 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x1 = in[1 * istride];
    sp_t x2 = in[6 * istride];
    sp_t x3 = in[11 * istride];
    sp_t x0 = in[16 * istride];

    a16 = sp_ntt_add(x0, x2, p);
    a17 = sp_ntt_sub(x0, x2, p);
    a18 = sp_ntt_add(x1, x3, p);
    a19 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a04;
    b4 = a08;
    b2 = a12;
    b3 = a16;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a00 = b0;
    a04 = b4;
    a08 = b3;
    a12 = b1;
    a16 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a05;
    b4 = a09;
    b2 = a13;
    b3 = a17;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a01 = b0;
    a05 = b4;
    a09 = b3;
    a13 = b1;
    a17 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a06;
    b4 = a10;
    b2 = a14;
    b3 = a18;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a02 = b0;
    a06 = b4;
    a10 = b3;
    a14 = b1;
    a18 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a07;
    b4 = a11;
    b2 = a15;
    b3 = a19;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b0 = sp_ntt_mul(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a03 = b0;
    a07 = b4;
    a11 = b3;
    a15 = b1;
    a19 = b2;
  }
  {
    sp_t p0 = sp_ntt_add(a00, a02, p);
    sp_t p1 = sp_ntt_sub(a00, a02, p);
    sp_t p2 = sp_ntt_add(a01, a03, p);
    sp_t p3 = sp_ntt_sub(a01, a03, p);

    out[0 * ostride] = p0;
    out[5 * ostride] = p2;
    out[10 * ostride] = p1;
    out[15 * ostride] = p3;
  }
  {
    sp_t p0 = sp_ntt_add(a04, a06, p);
    sp_t p1 = sp_ntt_sub(a04, a06, p);
    sp_t p2 = sp_ntt_add(a05, a07, p);
    sp_t p3 = sp_ntt_sub(a05, a07, p);

    out[1 * ostride] = p2;
    out[6 * ostride] = p1;
    out[11 * ostride] = p3;
    out[16 * ostride] = p0;
  }
  {
    sp_t p0 = sp_ntt_add(a08, a10, p);
    sp_t p1 = sp_ntt_sub(a08, a10, p);
    sp_t p2 = sp_ntt_add(a09, a11, p);
    sp_t p3 = sp_ntt_sub(a09, a11, p);

    out[2 * ostride] = p1;
    out[7 * ostride] = p3;
    out[12 * ostride] = p0;
    out[17 * ostride] = p2;
  }
  {
    sp_t p0 = sp_ntt_add(a12, a14, p);
    sp_t p1 = sp_ntt_sub(a12, a14, p);
    sp_t p2 = sp_ntt_add(a13, a15, p);
    sp_t p3 = sp_ntt_sub(a13, a15, p);

    out[3 * ostride] = p3;
    out[8 * ostride] = p0;
    out[13 * ostride] = p2;
    out[18 * ostride] = p1;
  }
  {
    sp_t p0 = sp_ntt_add(a16, a18, p);
    sp_t p1 = sp_ntt_sub(a16, a18, p);
    sp_t p2 = sp_ntt_add(a17, a19, p);
    sp_t p3 = sp_ntt_sub(a17, a19, p);

    out[4 * ostride] = p0;
    out[9 * ostride] = p2;
    out[14 * ostride] = p1;
    out[19 * ostride] = p3;
  }
}

static void
ntt20_twiddle_run_core(spv_t in, spv_size_t istride,
                spv_t out, spv_size_t ostride,
                spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, a04, 
       a05, a06, a07, a08, a09, 
       a10, a11, a12, a13, a14, 
       a15, a16, a17, a18, a19;

  {
    sp_t x0 = in[0 * istride];
    sp_t x1 = in[5 * istride];
    sp_t x2 = in[10 * istride];
    sp_t x3 = in[15 * istride];

    a00 = sp_ntt_add(x0, x2, p);
    a01 = sp_ntt_sub(x0, x2, p);
    a02 = sp_ntt_add(x1, x3, p);
    a03 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x0 = in[4 * istride];
    sp_t x1 = in[9 * istride];
    sp_t x2 = in[14 * istride];
    sp_t x3 = in[19 * istride];

    a04 = sp_ntt_add(x0, x2, p);
    a05 = sp_ntt_sub(x0, x2, p);
    a06 = sp_ntt_add(x1, x3, p);
    a07 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x3 = in[3 * istride];
    sp_t x0 = in[8 * istride];
    sp_t x1 = in[13 * istride];
    sp_t x2 = in[18 * istride];

    a08 = sp_ntt_add(x0, x2, p);
    a09 = sp_ntt_sub(x0, x2, p);
    a10 = sp_ntt_add(x1, x3, p);
    a11 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x2 = in[2 * istride];
    sp_t x3 = in[7 * istride];
    sp_t x0 = in[12 * istride];
    sp_t x1 = in[17 * istride];

    a12 = sp_ntt_add(x0, x2, p);
    a13 = sp_ntt_sub(x0, x2, p);
    a14 = sp_ntt_add(x1, x3, p);
    a15 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x1 = in[1 * istride];
    sp_t x2 = in[6 * istride];
    sp_t x3 = in[11 * istride];
    sp_t x0 = in[16 * istride];

    a16 = sp_ntt_add(x0, x2, p);
    a17 = sp_ntt_sub(x0, x2, p);
    a18 = sp_ntt_add(x1, x3, p);
    a19 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a04;
    b4 = a08;
    b2 = a12;
    b3 = a16;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a00 = b0;
    a04 = b4;
    a08 = b3;
    a12 = b1;
    a16 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a05;
    b4 = a09;
    b2 = a13;
    b3 = a17;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a01 = b0;
    a05 = b4;
    a09 = b3;
    a13 = b1;
    a17 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a06;
    b4 = a10;
    b2 = a14;
    b3 = a18;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a02 = b0;
    a06 = b4;
    a10 = b3;
    a14 = b1;
    a18 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a07;
    b4 = a11;
    b2 = a15;
    b3 = a19;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b0 = sp_ntt_mul(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a03 = b0;
    a07 = b4;
    a11 = b3;
    a15 = b1;
    a19 = b2;
  }
  {
    sp_t p0 = sp_ntt_add(a00, a02, p);
    sp_t p1 = sp_ntt_sub_partial(a00, a02, p);
    sp_t p2 = sp_ntt_add_partial(a01, a03, p);
    sp_t p3 = sp_ntt_sub_partial(a01, a03, p);

    p2 = sp_ntt_mul(p2, w[8], w[9], p);
    p1 = sp_ntt_mul(p1, w[18], w[19], p);
    p3 = sp_ntt_mul(p3, w[28], w[29], p);

    out[0 * ostride] = p0;
    out[5 * ostride] = p2;
    out[10 * ostride] = p1;
    out[15 * ostride] = p3;
  }
  {
    sp_t p0 = sp_ntt_add_partial(a04, a06, p);
    sp_t p1 = sp_ntt_sub_partial(a04, a06, p);
    sp_t p2 = sp_ntt_add_partial(a05, a07, p);
    sp_t p3 = sp_ntt_sub_partial(a05, a07, p);

    p2 = sp_ntt_mul(p2, w[0], w[1], p);
    p1 = sp_ntt_mul(p1, w[10], w[11], p);
    p3 = sp_ntt_mul(p3, w[20], w[21], p);
    p0 = sp_ntt_mul(p0, w[30], w[31], p);

    out[1 * ostride] = p2;
    out[6 * ostride] = p1;
    out[11 * ostride] = p3;
    out[16 * ostride] = p0;
  }
  {
    sp_t p0 = sp_ntt_add_partial(a08, a10, p);
    sp_t p1 = sp_ntt_sub_partial(a08, a10, p);
    sp_t p2 = sp_ntt_add_partial(a09, a11, p);
    sp_t p3 = sp_ntt_sub_partial(a09, a11, p);

    p1 = sp_ntt_mul(p1, w[2], w[3], p);
    p3 = sp_ntt_mul(p3, w[12], w[13], p);
    p0 = sp_ntt_mul(p0, w[22], w[23], p);
    p2 = sp_ntt_mul(p2, w[32], w[33], p);

    out[2 * ostride] = p1;
    out[7 * ostride] = p3;
    out[12 * ostride] = p0;
    out[17 * ostride] = p2;
  }
  {
    sp_t p0 = sp_ntt_add_partial(a12, a14, p);
    sp_t p1 = sp_ntt_sub_partial(a12, a14, p);
    sp_t p2 = sp_ntt_add_partial(a13, a15, p);
    sp_t p3 = sp_ntt_sub_partial(a13, a15, p);

    p3 = sp_ntt_mul(p3, w[4], w[5], p);
    p0 = sp_ntt_mul(p0, w[14], w[15], p);
    p2 = sp_ntt_mul(p2, w[24], w[25], p);
    p1 = sp_ntt_mul(p1, w[34], w[35], p);

    out[3 * ostride] = p3;
    out[8 * ostride] = p0;
    out[13 * ostride] = p2;
    out[18 * ostride] = p1;
  }
  {
    sp_t p0 = sp_ntt_add_partial(a16, a18, p);
    sp_t p1 = sp_ntt_sub_partial(a16, a18, p);
    sp_t p2 = sp_ntt_add_partial(a17, a19, p);
    sp_t p3 = sp_ntt_sub_partial(a17, a19, p);

    p0 = sp_ntt_mul(p0, w[6], w[7], p);
    p2 = sp_ntt_mul(p2, w[16], w[17], p);
    p1 = sp_ntt_mul(p1, w[26], w[27], p);
    p3 = sp_ntt_mul(p3, w[36], w[37], p);

    out[4 * ostride] = p0;
    out[9 * ostride] = p2;
    out[14 * ostride] = p1;
    out[19 * ostride] = p3;
  }
}

static void
ntt20_pfa_run_core(spv_t x, spv_size_t start,
          spv_size_t inc, spv_size_t n,
          sp_t p, spv_t ntt_const)
{
  spv_size_t j00, j01, j02, j03, j04, 
	     j05, j06, j07, j08, j09, 
	     j10, j11, j12, j13, j14, 
	     j15, j16, j17, j18, j19;

  sp_t a00, a01, a02, a03, a04, 
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
    sp_t x0 = x[j00];
    sp_t x1 = x[j05];
    sp_t x2 = x[j10];
    sp_t x3 = x[j15];

    a00 = sp_ntt_add(x0, x2, p);
    a01 = sp_ntt_sub(x0, x2, p);
    a02 = sp_ntt_add(x1, x3, p);
    a03 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x0 = x[j04];
    sp_t x1 = x[j09];
    sp_t x2 = x[j14];
    sp_t x3 = x[j19];

    a04 = sp_ntt_add(x0, x2, p);
    a05 = sp_ntt_sub(x0, x2, p);
    a06 = sp_ntt_add(x1, x3, p);
    a07 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x3 = x[j03];
    sp_t x0 = x[j08];
    sp_t x1 = x[j13];
    sp_t x2 = x[j18];

    a08 = sp_ntt_add(x0, x2, p);
    a09 = sp_ntt_sub(x0, x2, p);
    a10 = sp_ntt_add(x1, x3, p);
    a11 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x2 = x[j02];
    sp_t x3 = x[j07];
    sp_t x0 = x[j12];
    sp_t x1 = x[j17];

    a12 = sp_ntt_add(x0, x2, p);
    a13 = sp_ntt_sub(x0, x2, p);
    a14 = sp_ntt_add(x1, x3, p);
    a15 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t x1 = x[j01];
    sp_t x2 = x[j06];
    sp_t x3 = x[j11];
    sp_t x0 = x[j16];

    a16 = sp_ntt_add(x0, x2, p);
    a17 = sp_ntt_sub(x0, x2, p);
    a18 = sp_ntt_add(x1, x3, p);
    a19 = sp_ntt_sub(x1, x3, p);
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a00;
    b1 = a04;
    b4 = a08;
    b2 = a12;
    b3 = a16;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[1], ntt_const[NC+1], p);
    b2 = sp_ntt_mul(b2, ntt_const[2], ntt_const[NC+2], p);
    b3 = sp_ntt_mul(b3, ntt_const[3], ntt_const[NC+3], p);
    b4 = sp_ntt_mul(b4, ntt_const[4], ntt_const[NC+4], p);
    b5 = sp_ntt_mul(b5, ntt_const[5], ntt_const[NC+5], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a00 = b0;
    a04 = b4;
    a08 = b3;
    a12 = b1;
    a16 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a01;
    b1 = a05;
    b4 = a09;
    b2 = a13;
    b3 = a17;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[7], ntt_const[NC+7], p);
    b2 = sp_ntt_mul(b2, ntt_const[8], ntt_const[NC+8], p);
    b3 = sp_ntt_mul(b3, ntt_const[9], ntt_const[NC+9], p);
    b4 = sp_ntt_mul(b4, ntt_const[10], ntt_const[NC+10], p);
    b5 = sp_ntt_mul(b5, ntt_const[11], ntt_const[NC+11], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a01 = b0;
    a05 = b4;
    a09 = b3;
    a13 = b1;
    a17 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a02;
    b1 = a06;
    b4 = a10;
    b2 = a14;
    b3 = a18;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b1 = sp_ntt_mul(b1, ntt_const[13], ntt_const[NC+13], p);
    b2 = sp_ntt_mul(b2, ntt_const[14], ntt_const[NC+14], p);
    b3 = sp_ntt_mul(b3, ntt_const[15], ntt_const[NC+15], p);
    b4 = sp_ntt_mul(b4, ntt_const[16], ntt_const[NC+16], p);
    b5 = sp_ntt_mul(b5, ntt_const[17], ntt_const[NC+17], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a02 = b0;
    a06 = b4;
    a10 = b3;
    a14 = b1;
    a18 = b2;
  }
  {
    sp_t b0, b1, b2, b3, b4, b5;
    sp_t c1, c2, c3, c4;

    b0 = a03;
    b1 = a07;
    b4 = a11;
    b2 = a15;
    b3 = a19;

    c1 = sp_ntt_add(b1, b3, p);
    c3 = sp_ntt_sub(b1, b3, p);
    c2 = sp_ntt_add(b2, b4, p);
    c4 = sp_ntt_sub(b2, b4, p);

    b1 = sp_ntt_add(c1, c2, p);
    b2 = sp_ntt_sub_partial(c1, c2, p);
    b3 = c3;
    b4 = c4;
    b5 = sp_ntt_add_partial(c3, c4, p);

    b0 = sp_ntt_add(b0, b1, p);

    b0 = sp_ntt_mul(b0, ntt_const[18], ntt_const[NC+18], p);
    b1 = sp_ntt_mul(b1, ntt_const[19], ntt_const[NC+19], p);
    b2 = sp_ntt_mul(b2, ntt_const[20], ntt_const[NC+20], p);
    b3 = sp_ntt_mul(b3, ntt_const[21], ntt_const[NC+21], p);
    b4 = sp_ntt_mul(b4, ntt_const[22], ntt_const[NC+22], p);
    b5 = sp_ntt_mul(b5, ntt_const[23], ntt_const[NC+23], p);

    b1 = sp_ntt_add(b0, b1, p);

    c1 = sp_ntt_add(b1, b2, p);
    c2 = sp_ntt_sub(b1, b2, p);
    c3 = sp_ntt_add(b3, b5, p);
    c4 = sp_ntt_add(b4, b5, p);

    b1 = sp_ntt_add(c1, c3, p);
    b2 = sp_ntt_add(c2, c4, p);
    b3 = sp_ntt_sub(c1, c3, p);
    b4 = sp_ntt_sub(c2, c4, p);

    a03 = b0;
    a07 = b4;
    a11 = b3;
    a15 = b1;
    a19 = b2;
  }
  {
    sp_t p0 = sp_ntt_add(a00, a02, p);
    sp_t p1 = sp_ntt_sub(a00, a02, p);
    sp_t p2 = sp_ntt_add(a01, a03, p);
    sp_t p3 = sp_ntt_sub(a01, a03, p);

    x[j00] = p0;
    x[j05] = p2;
    x[j10] = p1;
    x[j15] = p3;
  }
  {
    sp_t p0 = sp_ntt_add(a04, a06, p);
    sp_t p1 = sp_ntt_sub(a04, a06, p);
    sp_t p2 = sp_ntt_add(a05, a07, p);
    sp_t p3 = sp_ntt_sub(a05, a07, p);

    x[j01] = p2;
    x[j06] = p1;
    x[j11] = p3;
    x[j16] = p0;
  }
  {
    sp_t p0 = sp_ntt_add(a08, a10, p);
    sp_t p1 = sp_ntt_sub(a08, a10, p);
    sp_t p2 = sp_ntt_add(a09, a11, p);
    sp_t p3 = sp_ntt_sub(a09, a11, p);

    x[j02] = p1;
    x[j07] = p3;
    x[j12] = p0;
    x[j17] = p2;
  }
  {
    sp_t p0 = sp_ntt_add(a12, a14, p);
    sp_t p1 = sp_ntt_sub(a12, a14, p);
    sp_t p2 = sp_ntt_add(a13, a15, p);
    sp_t p3 = sp_ntt_sub(a13, a15, p);

    x[j03] = p3;
    x[j08] = p0;
    x[j13] = p2;
    x[j18] = p1;
  }
  {
    sp_t p0 = sp_ntt_add(a16, a18, p);
    sp_t p1 = sp_ntt_sub(a16, a18, p);
    sp_t p2 = sp_ntt_add(a17, a19, p);
    sp_t p3 = sp_ntt_sub(a17, a19, p);

    x[j04] = p0;
    x[j09] = p2;
    x[j14] = p1;
    x[j19] = p3;
  }
}

DECLARE_CORE_ROUTINES(20)
