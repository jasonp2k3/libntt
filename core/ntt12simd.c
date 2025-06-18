#include "ntt/ntt-impl-simd.h"

#define NC 12

static const uint8_t ntt12_fixed_const[NC] = {1, 1, 1};

extern void X(ntt12_init)(spv_t out, sp_t p, sp_t d, sp_t primroot, 
                        sp_t order, sp_t perm);


static void
ntt12_run_core_simd(spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, 
            a04, a05, a06, a07, 
            a08, a09, a10, a11;

  {
    sp_simd_t x00 = sp_simd_gather(in + 0 * istride, idist, vsize);
    sp_simd_t x01 = sp_simd_gather(in + 4 * istride, idist, vsize);
    sp_simd_t x02 = sp_simd_gather(in + 8 * istride, idist, vsize);

    a01 = sp_ntt_add_simd(x01, x02, p);
    a02 = sp_ntt_sub_simd(x01, x02, p);

    a00 = sp_ntt_add_simd(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_gather(in + 3 * istride, idist, vsize);
    sp_simd_t x04 = sp_simd_gather(in + 7 * istride, idist, vsize);
    sp_simd_t x05 = sp_simd_gather(in +11 * istride, idist, vsize);

    a04 = sp_ntt_add_simd(x04, x05, p);
    a05 = sp_ntt_sub_simd(x04, x05, p);

    a03 = sp_ntt_add_simd(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_gather(in + 2 * istride, idist, vsize);
    sp_simd_t x06 = sp_simd_gather(in + 6 * istride, idist, vsize);
    sp_simd_t x07 = sp_simd_gather(in +10 * istride, idist, vsize);

    a07 = sp_ntt_add_simd(x07, x08, p);
    a08 = sp_ntt_sub_simd(x07, x08, p);

    a06 = sp_ntt_add_simd(x06, a07, p);
  }
  {
    sp_simd_t x10 = sp_simd_gather(in + 1 * istride, idist, vsize);
    sp_simd_t x11 = sp_simd_gather(in + 5 * istride, idist, vsize);
    sp_simd_t x09 = sp_simd_gather(in + 9 * istride, idist, vsize);

    a10 = sp_ntt_add_simd(x10, x11, p);
    a11 = sp_ntt_sub_simd(x10, x11, p);

    a09 = sp_ntt_add_simd(x09, a10, p);
  }
  {
    sp_simd_t x0 = a00;
    sp_simd_t x1 = a03;
    sp_simd_t x2 = a06;
    sp_simd_t x3 = a09;

    sp_simd_t t0 = sp_ntt_add_simd(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_simd(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_simd(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd(x1, x3, p);

    t3 = sp_ntt_mul_simd(t3, ntt_const[3], ntt_const[NC+3], p);

    x0 = sp_ntt_add_simd(t0, t1, p);
    x1 = sp_ntt_sub_simd(t0, t1, p);
    x2 = sp_ntt_add_simd(t2, t3, p);
    x3 = sp_ntt_sub_simd(t2, t3, p);

    a00 = x0;
    a03 = x2;
    a06 = x1;
    a09 = x3;
  }
  {
    sp_simd_t x0 = a01;
    sp_simd_t x1 = a04;
    sp_simd_t x2 = a07;
    sp_simd_t x3 = a10;

    sp_simd_t t0 = sp_ntt_add_partial_simd(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd(x1, x3, p);

    t0 = sp_ntt_mul_simd(t0, ntt_const[4], ntt_const[NC+4], p);
    t2 = sp_ntt_mul_simd(t2, ntt_const[5], ntt_const[NC+5], p);
    t1 = sp_ntt_mul_simd(t1, ntt_const[6], ntt_const[NC+6], p);
    t3 = sp_ntt_mul_simd(t3, ntt_const[7], ntt_const[NC+7], p);

    x0 = sp_ntt_add_simd(t0, t1, p);
    x1 = sp_ntt_sub_simd(t0, t1, p);
    x2 = sp_ntt_add_simd(t2, t3, p);
    x3 = sp_ntt_sub_simd(t2, t3, p);

    a01 = x0;
    a04 = x2;
    a07 = x1;
    a10 = x3;
  }
  {
    sp_simd_t x0 = a02;
    sp_simd_t x1 = a05;
    sp_simd_t x2 = a08;
    sp_simd_t x3 = a11;

    sp_simd_t t0 = sp_ntt_add_partial_simd(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd(x1, x3, p);

    t0 = sp_ntt_mul_simd(t0, ntt_const[8], ntt_const[NC+8], p);
    t2 = sp_ntt_mul_simd(t2, ntt_const[9], ntt_const[NC+9], p);
    t1 = sp_ntt_mul_simd(t1, ntt_const[10], ntt_const[NC+10], p);
    t3 = sp_ntt_mul_simd(t3, ntt_const[11], ntt_const[NC+11], p);

    x0 = sp_ntt_add_simd(t0, t1, p);
    x1 = sp_ntt_sub_simd(t0, t1, p);
    x2 = sp_ntt_add_simd(t2, t3, p);
    x3 = sp_ntt_sub_simd(t2, t3, p);

    a02 = x0;
    a05 = x2;
    a08 = x1;
    a11 = x3;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_simd(a01, a02, p);
    x02 = sp_ntt_sub_simd(a01, a02, p);

    sp_simd_scatter(x00, out + 0 * ostride, odist, vsize);
    sp_simd_scatter(x01, out + 4 * ostride, odist, vsize);
    sp_simd_scatter(x02, out + 8 * ostride, odist, vsize);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd(a04, a05, p);
    x05 = sp_ntt_sub_simd(a04, a05, p);

    sp_simd_scatter(x04, out + 1 * ostride, odist, vsize);
    sp_simd_scatter(x05, out + 5 * ostride, odist, vsize);
    sp_simd_scatter(x03, out + 9 * ostride, odist, vsize);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd(a07, a08, p);
    x08 = sp_ntt_sub_simd(a07, a08, p);

    sp_simd_scatter(x08, out + 2 * ostride, odist, vsize);
    sp_simd_scatter(x06, out + 6 * ostride, odist, vsize);
    sp_simd_scatter(x07, out +10 * ostride, odist, vsize);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd(a10, a11, p);
    x11 = sp_ntt_sub_simd(a10, a11, p);

    sp_simd_scatter(x09, out + 3 * ostride, odist, vsize);
    sp_simd_scatter(x10, out + 7 * ostride, odist, vsize);
    sp_simd_scatter(x11, out +11 * ostride, odist, vsize);
  }
}

static void
ntt12_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, 
            a04, a05, a06, a07, 
            a08, a09, a10, a11;

  {
    sp_simd_t x00 = sp_simd_load(in + 0 * istride);
    sp_simd_t x01 = sp_simd_load(in + 4 * istride);
    sp_simd_t x02 = sp_simd_load(in + 8 * istride);

    a01 = sp_ntt_add_simd0(x01, x02, p);
    a02 = sp_ntt_sub_simd0(x01, x02, p);

    a00 = sp_ntt_add_simd0(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_load(in + 3 * istride);
    sp_simd_t x04 = sp_simd_load(in + 7 * istride);
    sp_simd_t x05 = sp_simd_load(in +11 * istride);

    a04 = sp_ntt_add_simd0(x04, x05, p);
    a05 = sp_ntt_sub_simd0(x04, x05, p);

    a03 = sp_ntt_add_simd0(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_load(in + 2 * istride);
    sp_simd_t x06 = sp_simd_load(in + 6 * istride);
    sp_simd_t x07 = sp_simd_load(in +10 * istride);

    a07 = sp_ntt_add_simd0(x07, x08, p);
    a08 = sp_ntt_sub_simd0(x07, x08, p);

    a06 = sp_ntt_add_simd0(x06, a07, p);
  }
  {
    sp_simd_t x10 = sp_simd_load(in + 1 * istride);
    sp_simd_t x11 = sp_simd_load(in + 5 * istride);
    sp_simd_t x09 = sp_simd_load(in + 9 * istride);

    a10 = sp_ntt_add_simd0(x10, x11, p);
    a11 = sp_ntt_sub_simd0(x10, x11, p);

    a09 = sp_ntt_add_simd0(x09, a10, p);
  }
  {
    sp_simd_t x0 = a00;
    sp_simd_t x1 = a03;
    sp_simd_t x2 = a06;
    sp_simd_t x3 = a09;

    sp_simd_t t0 = sp_ntt_add_simd0(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_simd0(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_simd0(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd0(x1, x3, p);

    t3 = sp_ntt_mul_simd0(t3, c + 2*3, p);

    x0 = sp_ntt_add_simd0(t0, t1, p);
    x1 = sp_ntt_sub_simd0(t0, t1, p);
    x2 = sp_ntt_add_simd0(t2, t3, p);
    x3 = sp_ntt_sub_simd0(t2, t3, p);

    a00 = x0;
    a03 = x2;
    a06 = x1;
    a09 = x3;
  }
  {
    sp_simd_t x0 = a01;
    sp_simd_t x1 = a04;
    sp_simd_t x2 = a07;
    sp_simd_t x3 = a10;

    sp_simd_t t0 = sp_ntt_add_partial_simd0(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd0(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd0(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd0(x1, x3, p);

    t0 = sp_ntt_mul_simd0(t0, c + 2*4, p);
    t2 = sp_ntt_mul_simd0(t2, c + 2*5, p);
    t1 = sp_ntt_mul_simd0(t1, c + 2*6, p);
    t3 = sp_ntt_mul_simd0(t3, c + 2*7, p);

    x0 = sp_ntt_add_simd0(t0, t1, p);
    x1 = sp_ntt_sub_simd0(t0, t1, p);
    x2 = sp_ntt_add_simd0(t2, t3, p);
    x3 = sp_ntt_sub_simd0(t2, t3, p);

    a01 = x0;
    a04 = x2;
    a07 = x1;
    a10 = x3;
  }
  {
    sp_simd_t x0 = a02;
    sp_simd_t x1 = a05;
    sp_simd_t x2 = a08;
    sp_simd_t x3 = a11;

    sp_simd_t t0 = sp_ntt_add_partial_simd0(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd0(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd0(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd0(x1, x3, p);

    t0 = sp_ntt_mul_simd0(t0, c + 2*8, p);
    t2 = sp_ntt_mul_simd0(t2, c + 2*9, p);
    t1 = sp_ntt_mul_simd0(t1, c + 2*10, p);
    t3 = sp_ntt_mul_simd0(t3, c + 2*11, p);

    x0 = sp_ntt_add_simd0(t0, t1, p);
    x1 = sp_ntt_sub_simd0(t0, t1, p);
    x2 = sp_ntt_add_simd0(t2, t3, p);
    x3 = sp_ntt_sub_simd0(t2, t3, p);

    a02 = x0;
    a05 = x2;
    a08 = x1;
    a11 = x3;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd0(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_simd0(a01, a02, p);
    x02 = sp_ntt_sub_simd0(a01, a02, p);

    sp_simd_store(x00, out + 0 * ostride);
    sp_simd_store(x01, out + 4 * ostride);
    sp_simd_store(x02, out + 8 * ostride);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd0(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd0(a04, a05, p);
    x05 = sp_ntt_sub_simd0(a04, a05, p);

    sp_simd_store(x04, out + 1 * ostride);
    sp_simd_store(x05, out + 5 * ostride);
    sp_simd_store(x03, out + 9 * ostride);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd0(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd0(a07, a08, p);
    x08 = sp_ntt_sub_simd0(a07, a08, p);

    sp_simd_store(x08, out + 2 * ostride);
    sp_simd_store(x06, out + 6 * ostride);
    sp_simd_store(x07, out +10 * ostride);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd0(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd0(a10, a11, p);
    x11 = sp_ntt_sub_simd0(a10, a11, p);

    sp_simd_store(x09, out + 3 * ostride);
    sp_simd_store(x10, out + 7 * ostride);
    sp_simd_store(x11, out +11 * ostride);
  }
}


static void
ntt12_twiddle_run_core_simd(
        spv_t in, spv_size_t istride, spv_size_t idist,
		spv_t out, spv_size_t ostride, spv_size_t odist,
		sp_simd_t *w, sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  sp_simd_t a00, a01, a02, a03, 
            a04, a05, a06, a07, 
            a08, a09, a10, a11;

  {
    sp_simd_t x00 = sp_simd_gather(in + 0 * istride, idist, vsize);
    sp_simd_t x01 = sp_simd_gather(in + 4 * istride, idist, vsize);
    sp_simd_t x02 = sp_simd_gather(in + 8 * istride, idist, vsize);

    a01 = sp_ntt_add_simd(x01, x02, p);
    a02 = sp_ntt_sub_simd(x01, x02, p);

    a00 = sp_ntt_add_simd(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_gather(in + 3 * istride, idist, vsize);
    sp_simd_t x04 = sp_simd_gather(in + 7 * istride, idist, vsize);
    sp_simd_t x05 = sp_simd_gather(in +11 * istride, idist, vsize);

    a04 = sp_ntt_add_simd(x04, x05, p);
    a05 = sp_ntt_sub_simd(x04, x05, p);

    a03 = sp_ntt_add_simd(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_gather(in + 2 * istride, idist, vsize);
    sp_simd_t x06 = sp_simd_gather(in + 6 * istride, idist, vsize);
    sp_simd_t x07 = sp_simd_gather(in +10 * istride, idist, vsize);

    a07 = sp_ntt_add_simd(x07, x08, p);
    a08 = sp_ntt_sub_simd(x07, x08, p);

    a06 = sp_ntt_add_simd(x06, a07, p);
  }
  {
    sp_simd_t x10 = sp_simd_gather(in + 1 * istride, idist, vsize);
    sp_simd_t x11 = sp_simd_gather(in + 5 * istride, idist, vsize);
    sp_simd_t x09 = sp_simd_gather(in + 9 * istride, idist, vsize);

    a10 = sp_ntt_add_simd(x10, x11, p);
    a11 = sp_ntt_sub_simd(x10, x11, p);

    a09 = sp_ntt_add_simd(x09, a10, p);
  }
  {
    sp_simd_t x0 = a00;
    sp_simd_t x1 = a03;
    sp_simd_t x2 = a06;
    sp_simd_t x3 = a09;

    sp_simd_t t0 = sp_ntt_add_simd(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_simd(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_simd(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd(x1, x3, p);

    t3 = sp_ntt_mul_simd(t3, ntt_const[3], ntt_const[NC+3], p);

    x0 = sp_ntt_add_simd(t0, t1, p);
    x1 = sp_ntt_sub_simd(t0, t1, p);
    x2 = sp_ntt_add_simd(t2, t3, p);
    x3 = sp_ntt_sub_simd(t2, t3, p);

    a00 = x0;
    a03 = x2;
    a06 = x1;
    a09 = x3;
  }
  {
    sp_simd_t x0 = a01;
    sp_simd_t x1 = a04;
    sp_simd_t x2 = a07;
    sp_simd_t x3 = a10;

    sp_simd_t t0 = sp_ntt_add_partial_simd(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd(x1, x3, p);

    t0 = sp_ntt_mul_simd(t0, ntt_const[4], ntt_const[NC+4], p);
    t2 = sp_ntt_mul_simd(t2, ntt_const[5], ntt_const[NC+5], p);
    t1 = sp_ntt_mul_simd(t1, ntt_const[6], ntt_const[NC+6], p);
    t3 = sp_ntt_mul_simd(t3, ntt_const[7], ntt_const[NC+7], p);

    x0 = sp_ntt_add_simd(t0, t1, p);
    x1 = sp_ntt_sub_simd(t0, t1, p);
    x2 = sp_ntt_add_simd(t2, t3, p);
    x3 = sp_ntt_sub_simd(t2, t3, p);

    a01 = x0;
    a04 = x2;
    a07 = x1;
    a10 = x3;
  }
  {
    sp_simd_t x0 = a02;
    sp_simd_t x1 = a05;
    sp_simd_t x2 = a08;
    sp_simd_t x3 = a11;

    sp_simd_t t0 = sp_ntt_add_partial_simd(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd(x1, x3, p);

    t0 = sp_ntt_mul_simd(t0, ntt_const[8], ntt_const[NC+8], p);
    t2 = sp_ntt_mul_simd(t2, ntt_const[9], ntt_const[NC+9], p);
    t1 = sp_ntt_mul_simd(t1, ntt_const[10], ntt_const[NC+10], p);
    t3 = sp_ntt_mul_simd(t3, ntt_const[11], ntt_const[NC+11], p);

    x0 = sp_ntt_add_simd(t0, t1, p);
    x1 = sp_ntt_sub_simd(t0, t1, p);
    x2 = sp_ntt_add_simd(t2, t3, p);
    x3 = sp_ntt_sub_simd(t2, t3, p);

    a02 = x0;
    a05 = x2;
    a08 = x1;
    a11 = x3;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_partial_simd(a01, a02, p);
    x02 = sp_ntt_sub_partial_simd(a01, a02, p);

    x01 = sp_ntt_twiddle_mul_simd(x01, w + 6, p);
    x02 = sp_ntt_twiddle_mul_simd(x02, w + 14, p);

    sp_simd_scatter(x00, out + 0 * ostride, odist, vsize);
    sp_simd_scatter(x01, out + 4 * ostride, odist, vsize);
    sp_simd_scatter(x02, out + 8 * ostride, odist, vsize);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_partial_simd(a04, a05, p);
    x05 = sp_ntt_sub_partial_simd(a04, a05, p);

    x04 = sp_ntt_twiddle_mul_simd(x04, w + 0, p);
    x05 = sp_ntt_twiddle_mul_simd(x05, w + 8, p);
    x03 = sp_ntt_twiddle_mul_simd(x03, w + 16, p);

    sp_simd_scatter(x04, out + 1 * ostride, odist, vsize);
    sp_simd_scatter(x05, out + 5 * ostride, odist, vsize);
    sp_simd_scatter(x03, out + 9 * ostride, odist, vsize);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_partial_simd(a07, a08, p);
    x08 = sp_ntt_sub_partial_simd(a07, a08, p);

    x08 = sp_ntt_twiddle_mul_simd(x08, w + 2, p);
    x06 = sp_ntt_twiddle_mul_simd(x06, w + 10, p);
    x07 = sp_ntt_twiddle_mul_simd(x07, w + 18, p);

    sp_simd_scatter(x08, out + 2 * ostride, odist, vsize);
    sp_simd_scatter(x06, out + 6 * ostride, odist, vsize);
    sp_simd_scatter(x07, out +10 * ostride, odist, vsize);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_partial_simd(a10, a11, p);
    x11 = sp_ntt_sub_partial_simd(a10, a11, p);

    x09 = sp_ntt_twiddle_mul_simd(x09, w + 4, p);
    x10 = sp_ntt_twiddle_mul_simd(x10, w + 12, p);
    x11 = sp_ntt_twiddle_mul_simd(x11, w + 20, p);

    sp_simd_scatter(x09, out + 3 * ostride, odist, vsize);
    sp_simd_scatter(x10, out + 7 * ostride, odist, vsize);
    sp_simd_scatter(x11, out +11 * ostride, odist, vsize);
  }
}

static void
ntt12_twiddle_run_core_simd_interleaved(
		spv_t in, spv_size_t istride,
		spv_t out, spv_size_t ostride,
		sp_simd_t *w, sp_simd_t p, sp_simd_t * c)
{
  sp_simd_t a00, a01, a02, a03, 
            a04, a05, a06, a07, 
            a08, a09, a10, a11;

  {
    sp_simd_t x00 = sp_simd_load(in + 0 * istride);
    sp_simd_t x01 = sp_simd_load(in + 4 * istride);
    sp_simd_t x02 = sp_simd_load(in + 8 * istride);

    a01 = sp_ntt_add_simd0(x01, x02, p);
    a02 = sp_ntt_sub_simd0(x01, x02, p);

    a00 = sp_ntt_add_simd0(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_load(in + 3 * istride);
    sp_simd_t x04 = sp_simd_load(in + 7 * istride);
    sp_simd_t x05 = sp_simd_load(in +11 * istride);

    a04 = sp_ntt_add_simd0(x04, x05, p);
    a05 = sp_ntt_sub_simd0(x04, x05, p);

    a03 = sp_ntt_add_simd0(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_load(in + 2 * istride);
    sp_simd_t x06 = sp_simd_load(in + 6 * istride);
    sp_simd_t x07 = sp_simd_load(in +10 * istride);

    a07 = sp_ntt_add_simd0(x07, x08, p);
    a08 = sp_ntt_sub_simd0(x07, x08, p);

    a06 = sp_ntt_add_simd0(x06, a07, p);
  }
  {
    sp_simd_t x10 = sp_simd_load(in + 1 * istride);
    sp_simd_t x11 = sp_simd_load(in + 5 * istride);
    sp_simd_t x09 = sp_simd_load(in + 9 * istride);

    a10 = sp_ntt_add_simd0(x10, x11, p);
    a11 = sp_ntt_sub_simd0(x10, x11, p);

    a09 = sp_ntt_add_simd0(x09, a10, p);
  }
  {
    sp_simd_t x0 = a00;
    sp_simd_t x1 = a03;
    sp_simd_t x2 = a06;
    sp_simd_t x3 = a09;

    sp_simd_t t0 = sp_ntt_add_simd0(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_simd0(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_simd0(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd0(x1, x3, p);

    t3 = sp_ntt_mul_simd0(t3, c + 2*3, p);

    x0 = sp_ntt_add_simd0(t0, t1, p);
    x1 = sp_ntt_sub_simd0(t0, t1, p);
    x2 = sp_ntt_add_simd0(t2, t3, p);
    x3 = sp_ntt_sub_simd0(t2, t3, p);

    a00 = x0;
    a03 = x2;
    a06 = x1;
    a09 = x3;
  }
  {
    sp_simd_t x0 = a01;
    sp_simd_t x1 = a04;
    sp_simd_t x2 = a07;
    sp_simd_t x3 = a10;

    sp_simd_t t0 = sp_ntt_add_partial_simd0(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd0(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd0(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd0(x1, x3, p);

    t0 = sp_ntt_mul_simd0(t0, c + 2*4, p);
    t2 = sp_ntt_mul_simd0(t2, c + 2*5, p);
    t1 = sp_ntt_mul_simd0(t1, c + 2*6, p);
    t3 = sp_ntt_mul_simd0(t3, c + 2*7, p);

    x0 = sp_ntt_add_simd0(t0, t1, p);
    x1 = sp_ntt_sub_simd0(t0, t1, p);
    x2 = sp_ntt_add_simd0(t2, t3, p);
    x3 = sp_ntt_sub_simd0(t2, t3, p);

    a01 = x0;
    a04 = x2;
    a07 = x1;
    a10 = x3;
  }
  {
    sp_simd_t x0 = a02;
    sp_simd_t x1 = a05;
    sp_simd_t x2 = a08;
    sp_simd_t x3 = a11;

    sp_simd_t t0 = sp_ntt_add_partial_simd0(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd0(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd0(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd0(x1, x3, p);

    t0 = sp_ntt_mul_simd0(t0, c + 2*8, p);
    t2 = sp_ntt_mul_simd0(t2, c + 2*9, p);
    t1 = sp_ntt_mul_simd0(t1, c + 2*10, p);
    t3 = sp_ntt_mul_simd0(t3, c + 2*11, p);

    x0 = sp_ntt_add_simd0(t0, t1, p);
    x1 = sp_ntt_sub_simd0(t0, t1, p);
    x2 = sp_ntt_add_simd0(t2, t3, p);
    x3 = sp_ntt_sub_simd0(t2, t3, p);

    a02 = x0;
    a05 = x2;
    a08 = x1;
    a11 = x3;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd0(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_partial_simd0(a01, a02, p);
    x02 = sp_ntt_sub_partial_simd0(a01, a02, p);

    x01 = sp_ntt_twiddle_mul_simd0(x01, w + 6, p);
    x02 = sp_ntt_twiddle_mul_simd0(x02, w + 14, p);

    sp_simd_store(x00, out + 0 * ostride);
    sp_simd_store(x01, out + 4 * ostride);
    sp_simd_store(x02, out + 8 * ostride);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd0(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_partial_simd0(a04, a05, p);
    x05 = sp_ntt_sub_partial_simd0(a04, a05, p);

    x04 = sp_ntt_twiddle_mul_simd0(x04, w + 0, p);
    x05 = sp_ntt_twiddle_mul_simd0(x05, w + 8, p);
    x03 = sp_ntt_twiddle_mul_simd0(x03, w + 16, p);

    sp_simd_store(x04, out + 1 * ostride);
    sp_simd_store(x05, out + 5 * ostride);
    sp_simd_store(x03, out + 9 * ostride);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd0(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_partial_simd0(a07, a08, p);
    x08 = sp_ntt_sub_partial_simd0(a07, a08, p);

    x08 = sp_ntt_twiddle_mul_simd0(x08, w + 2, p);
    x06 = sp_ntt_twiddle_mul_simd0(x06, w + 10, p);
    x07 = sp_ntt_twiddle_mul_simd0(x07, w + 18, p);

    sp_simd_store(x08, out + 2 * ostride);
    sp_simd_store(x06, out + 6 * ostride);
    sp_simd_store(x07, out +10 * ostride);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd0(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_partial_simd0(a10, a11, p);
    x11 = sp_ntt_sub_partial_simd0(a10, a11, p);

    x09 = sp_ntt_twiddle_mul_simd0(x09, w + 4, p);
    x10 = sp_ntt_twiddle_mul_simd0(x10, w + 12, p);
    x11 = sp_ntt_twiddle_mul_simd0(x11, w + 20, p);

    sp_simd_store(x09, out + 3 * ostride);
    sp_simd_store(x10, out + 7 * ostride);
    sp_simd_store(x11, out +11 * ostride);
  }
}


static void
ntt12_pfa_run_core_simd(spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t inc2, spv_size_t n,
	  sp_t p, spv_t ntt_const, spv_size_t vsize)
{
  spv_size_t j00, j01, j02, j03, 
             j04, j05, j06, j07, 
             j08, j09, j10, j11;

  sp_simd_t a00, a01, a02, a03, 
            a04, a05, a06, a07, 
            a08, a09, a10, a11;

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

  {
    sp_simd_t x00 = sp_simd_pfa_gather(x, j00, inc2, n, vsize);
    sp_simd_t x01 = sp_simd_pfa_gather(x, j04, inc2, n, vsize);
    sp_simd_t x02 = sp_simd_pfa_gather(x, j08, inc2, n, vsize);

    a01 = sp_ntt_add_simd(x01, x02, p);
    a02 = sp_ntt_sub_simd(x01, x02, p);

    a00 = sp_ntt_add_simd(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_pfa_gather(x, j03, inc2, n, vsize);
    sp_simd_t x04 = sp_simd_pfa_gather(x, j07, inc2, n, vsize);
    sp_simd_t x05 = sp_simd_pfa_gather(x, j11, inc2, n, vsize);

    a04 = sp_ntt_add_simd(x04, x05, p);
    a05 = sp_ntt_sub_simd(x04, x05, p);

    a03 = sp_ntt_add_simd(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_pfa_gather(x, j02, inc2, n, vsize);
    sp_simd_t x06 = sp_simd_pfa_gather(x, j06, inc2, n, vsize);
    sp_simd_t x07 = sp_simd_pfa_gather(x, j10, inc2, n, vsize);

    a07 = sp_ntt_add_simd(x07, x08, p);
    a08 = sp_ntt_sub_simd(x07, x08, p);

    a06 = sp_ntt_add_simd(x06, a07, p);
  }
  {
    sp_simd_t x10 = sp_simd_pfa_gather(x, j01, inc2, n, vsize);
    sp_simd_t x11 = sp_simd_pfa_gather(x, j05, inc2, n, vsize);
    sp_simd_t x09 = sp_simd_pfa_gather(x, j09, inc2, n, vsize);

    a10 = sp_ntt_add_simd(x10, x11, p);
    a11 = sp_ntt_sub_simd(x10, x11, p);

    a09 = sp_ntt_add_simd(x09, a10, p);
  }
  {
    sp_simd_t x0 = a00;
    sp_simd_t x1 = a03;
    sp_simd_t x2 = a06;
    sp_simd_t x3 = a09;

    sp_simd_t t0 = sp_ntt_add_simd(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_simd(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_simd(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd(x1, x3, p);

    t3 = sp_ntt_mul_simd(t3, ntt_const[3], ntt_const[NC+3], p);

    x0 = sp_ntt_add_simd(t0, t1, p);
    x1 = sp_ntt_sub_simd(t0, t1, p);
    x2 = sp_ntt_add_simd(t2, t3, p);
    x3 = sp_ntt_sub_simd(t2, t3, p);

    a00 = x0;
    a03 = x2;
    a06 = x1;
    a09 = x3;
  }
  {
    sp_simd_t x0 = a01;
    sp_simd_t x1 = a04;
    sp_simd_t x2 = a07;
    sp_simd_t x3 = a10;

    sp_simd_t t0 = sp_ntt_add_partial_simd(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd(x1, x3, p);

    t0 = sp_ntt_mul_simd(t0, ntt_const[4], ntt_const[NC+4], p);
    t2 = sp_ntt_mul_simd(t2, ntt_const[5], ntt_const[NC+5], p);
    t1 = sp_ntt_mul_simd(t1, ntt_const[6], ntt_const[NC+6], p);
    t3 = sp_ntt_mul_simd(t3, ntt_const[7], ntt_const[NC+7], p);

    x0 = sp_ntt_add_simd(t0, t1, p);
    x1 = sp_ntt_sub_simd(t0, t1, p);
    x2 = sp_ntt_add_simd(t2, t3, p);
    x3 = sp_ntt_sub_simd(t2, t3, p);

    a01 = x0;
    a04 = x2;
    a07 = x1;
    a10 = x3;
  }
  {
    sp_simd_t x0 = a02;
    sp_simd_t x1 = a05;
    sp_simd_t x2 = a08;
    sp_simd_t x3 = a11;

    sp_simd_t t0 = sp_ntt_add_partial_simd(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd(x1, x3, p);

    t0 = sp_ntt_mul_simd(t0, ntt_const[8], ntt_const[NC+8], p);
    t2 = sp_ntt_mul_simd(t2, ntt_const[9], ntt_const[NC+9], p);
    t1 = sp_ntt_mul_simd(t1, ntt_const[10], ntt_const[NC+10], p);
    t3 = sp_ntt_mul_simd(t3, ntt_const[11], ntt_const[NC+11], p);

    x0 = sp_ntt_add_simd(t0, t1, p);
    x1 = sp_ntt_sub_simd(t0, t1, p);
    x2 = sp_ntt_add_simd(t2, t3, p);
    x3 = sp_ntt_sub_simd(t2, t3, p);

    a02 = x0;
    a05 = x2;
    a08 = x1;
    a11 = x3;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_simd(a01, a02, p);
    x02 = sp_ntt_sub_simd(a01, a02, p);

    sp_simd_pfa_scatter(x00, x, j00, inc2, n, vsize);
    sp_simd_pfa_scatter(x01, x, j04, inc2, n, vsize);
    sp_simd_pfa_scatter(x02, x, j08, inc2, n, vsize);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd(a04, a05, p);
    x05 = sp_ntt_sub_simd(a04, a05, p);

    sp_simd_pfa_scatter(x04, x, j01, inc2, n, vsize);
    sp_simd_pfa_scatter(x05, x, j05, inc2, n, vsize);
    sp_simd_pfa_scatter(x03, x, j09, inc2, n, vsize);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd(a07, a08, p);
    x08 = sp_ntt_sub_simd(a07, a08, p);

    sp_simd_pfa_scatter(x08, x, j02, inc2, n, vsize);
    sp_simd_pfa_scatter(x06, x, j06, inc2, n, vsize);
    sp_simd_pfa_scatter(x07, x, j10, inc2, n, vsize);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd(a10, a11, p);
    x11 = sp_ntt_sub_simd(a10, a11, p);

    sp_simd_pfa_scatter(x09, x, j03, inc2, n, vsize);
    sp_simd_pfa_scatter(x10, x, j07, inc2, n, vsize);
    sp_simd_pfa_scatter(x11, x, j11, inc2, n, vsize);
  }
}

static void
ntt12_pfa_run_core_simd_interleaved(
	  spv_t x, spv_size_t start,
	  spv_size_t inc, spv_size_t n,
	  sp_simd_t p, sp_simd_t * c)
{
  spv_size_t j00, j01, j02, j03, 
             j04, j05, j06, j07, 
             j08, j09, j10, j11;

  sp_simd_t a00, a01, a02, a03, 
            a04, a05, a06, a07, 
            a08, a09, a10, a11;

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

  {
    sp_simd_t x00 = sp_simd_load(x + j00);
    sp_simd_t x01 = sp_simd_load(x + j04);
    sp_simd_t x02 = sp_simd_load(x + j08);

    a01 = sp_ntt_add_simd0(x01, x02, p);
    a02 = sp_ntt_sub_simd0(x01, x02, p);

    a00 = sp_ntt_add_simd0(x00, a01, p);
  }
  {
    sp_simd_t x03 = sp_simd_load(x + j03);
    sp_simd_t x04 = sp_simd_load(x + j07);
    sp_simd_t x05 = sp_simd_load(x + j11);

    a04 = sp_ntt_add_simd0(x04, x05, p);
    a05 = sp_ntt_sub_simd0(x04, x05, p);

    a03 = sp_ntt_add_simd0(x03, a04, p);
  }
  {
    sp_simd_t x08 = sp_simd_load(x + j02);
    sp_simd_t x06 = sp_simd_load(x + j06);
    sp_simd_t x07 = sp_simd_load(x + j10);

    a07 = sp_ntt_add_simd0(x07, x08, p);
    a08 = sp_ntt_sub_simd0(x07, x08, p);

    a06 = sp_ntt_add_simd0(x06, a07, p);
  }
  {
    sp_simd_t x10 = sp_simd_load(x + j01);
    sp_simd_t x11 = sp_simd_load(x + j05);
    sp_simd_t x09 = sp_simd_load(x + j09);

    a10 = sp_ntt_add_simd0(x10, x11, p);
    a11 = sp_ntt_sub_simd0(x10, x11, p);

    a09 = sp_ntt_add_simd0(x09, a10, p);
  }
  {
    sp_simd_t x0 = a00;
    sp_simd_t x1 = a03;
    sp_simd_t x2 = a06;
    sp_simd_t x3 = a09;

    sp_simd_t t0 = sp_ntt_add_simd0(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_simd0(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_simd0(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd0(x1, x3, p);

    t3 = sp_ntt_mul_simd0(t3, c + 2*3, p);

    x0 = sp_ntt_add_simd0(t0, t1, p);
    x1 = sp_ntt_sub_simd0(t0, t1, p);
    x2 = sp_ntt_add_simd0(t2, t3, p);
    x3 = sp_ntt_sub_simd0(t2, t3, p);

    a00 = x0;
    a03 = x2;
    a06 = x1;
    a09 = x3;
  }
  {
    sp_simd_t x0 = a01;
    sp_simd_t x1 = a04;
    sp_simd_t x2 = a07;
    sp_simd_t x3 = a10;

    sp_simd_t t0 = sp_ntt_add_partial_simd0(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd0(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd0(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd0(x1, x3, p);

    t0 = sp_ntt_mul_simd0(t0, c + 2*4, p);
    t2 = sp_ntt_mul_simd0(t2, c + 2*5, p);
    t1 = sp_ntt_mul_simd0(t1, c + 2*6, p);
    t3 = sp_ntt_mul_simd0(t3, c + 2*7, p);

    x0 = sp_ntt_add_simd0(t0, t1, p);
    x1 = sp_ntt_sub_simd0(t0, t1, p);
    x2 = sp_ntt_add_simd0(t2, t3, p);
    x3 = sp_ntt_sub_simd0(t2, t3, p);

    a01 = x0;
    a04 = x2;
    a07 = x1;
    a10 = x3;
  }
  {
    sp_simd_t x0 = a02;
    sp_simd_t x1 = a05;
    sp_simd_t x2 = a08;
    sp_simd_t x3 = a11;

    sp_simd_t t0 = sp_ntt_add_partial_simd0(x0, x2, p);
    sp_simd_t t2 = sp_ntt_sub_partial_simd0(x0, x2, p);
    sp_simd_t t1 = sp_ntt_add_partial_simd0(x1, x3, p);
    sp_simd_t t3 = sp_ntt_sub_partial_simd0(x1, x3, p);

    t0 = sp_ntt_mul_simd0(t0, c + 2*8, p);
    t2 = sp_ntt_mul_simd0(t2, c + 2*9, p);
    t1 = sp_ntt_mul_simd0(t1, c + 2*10, p);
    t3 = sp_ntt_mul_simd0(t3, c + 2*11, p);

    x0 = sp_ntt_add_simd0(t0, t1, p);
    x1 = sp_ntt_sub_simd0(t0, t1, p);
    x2 = sp_ntt_add_simd0(t2, t3, p);
    x3 = sp_ntt_sub_simd0(t2, t3, p);

    a02 = x0;
    a05 = x2;
    a08 = x1;
    a11 = x3;
  }
  {
    sp_simd_t x00, x01, x02;

    a01 = sp_ntt_add_simd0(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_simd0(a01, a02, p);
    x02 = sp_ntt_sub_simd0(a01, a02, p);

    sp_simd_store(x00, x + j00);
    sp_simd_store(x01, x + j04);
    sp_simd_store(x02, x + j08);
  }
  {
    sp_simd_t x03, x04, x05;

    a04 = sp_ntt_add_simd0(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_simd0(a04, a05, p);
    x05 = sp_ntt_sub_simd0(a04, a05, p);

    sp_simd_store(x04, x + j01);
    sp_simd_store(x05, x + j05);
    sp_simd_store(x03, x + j09);
  }
  {
    sp_simd_t x06, x07, x08;

    a07 = sp_ntt_add_simd0(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_simd0(a07, a08, p);
    x08 = sp_ntt_sub_simd0(a07, a08, p);

    sp_simd_store(x08, x + j02);
    sp_simd_store(x06, x + j06);
    sp_simd_store(x07, x + j10);
  }
  {
    sp_simd_t x09, x10, x11;

    a10 = sp_ntt_add_simd0(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_simd0(a10, a11, p);
    x11 = sp_ntt_sub_simd0(a10, a11, p);

    sp_simd_store(x09, x + j03);
    sp_simd_store(x10, x + j07);
    sp_simd_store(x11, x + j11);
  }
}

DECLARE_CORE_ROUTINES_SIMD(12)
