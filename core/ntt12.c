#include "ntt/ntt-impl-scalar.h"

#define NC 12

static const uint8_t ntt12_fixed_const[NC] = {1, 1, 1};

void
X(ntt12_init)(spv_t out, sp_t p, sp_t d,
          sp_t primroot, sp_t order, sp_t perm)
{
  X(nttdata_init_generic)(&X(ntt12_config), out, p, d, primroot, order, perm);
}

static void 
ntt12_run_core(spv_t in, spv_size_t istride,
                spv_t out, spv_size_t ostride,
                sp_t p, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, 
       a04, a05, a06, a07,
       a08, a09, a10, a11;

  {
    sp_t x00 = in[ 0 * istride];
    sp_t x01 = in[ 4 * istride];
    sp_t x02 = in[ 8 * istride];

    a01 = sp_ntt_add(x01, x02, p);
    a02 = sp_ntt_sub(x01, x02, p);

    a00 = sp_ntt_add(x00, a01, p);
  }
  {
    sp_t x03 = in[ 3 * istride];
    sp_t x04 = in[ 7 * istride];
    sp_t x05 = in[11 * istride];

    a04 = sp_ntt_add(x04, x05, p);
    a05 = sp_ntt_sub(x04, x05, p);

    a03 = sp_ntt_add(x03, a04, p);
  }
  {
    sp_t x08 = in[ 2 * istride];
    sp_t x06 = in[ 6 * istride];
    sp_t x07 = in[10 * istride];

    a07 = sp_ntt_add(x07, x08, p);
    a08 = sp_ntt_sub(x07, x08, p);

    a06 = sp_ntt_add(x06, a07, p);
  }
  {
    sp_t x10 = in[ 1 * istride];
    sp_t x11 = in[ 5 * istride];
    sp_t x09 = in[ 9 * istride];

    a10 = sp_ntt_add(x10, x11, p);
    a11 = sp_ntt_sub(x10, x11, p);

    a09 = sp_ntt_add(x09, a10, p);
  }
  {
    sp_t x0 = a00;
    sp_t x1 = a03;
    sp_t x2 = a06;
    sp_t x3 = a09;

    sp_t t0 = sp_ntt_add(x0, x2, p);
    sp_t t2 = sp_ntt_sub(x0, x2, p);
    sp_t t1 = sp_ntt_add(x1, x3, p);
    sp_t t3 = sp_ntt_sub_partial(x1, x3, p);

    t3 = sp_ntt_mul(t3, ntt_const[3], ntt_const[NC+3], p);

    x0 = sp_ntt_add(t0, t1, p);
    x1 = sp_ntt_sub(t0, t1, p);
    x2 = sp_ntt_add(t2, t3, p);
    x3 = sp_ntt_sub(t2, t3, p);

    a00 = x0;
    a03 = x2;
    a06 = x1;
    a09 = x3;
  }
  {
    sp_t x0 = a01;
    sp_t x1 = a04;
    sp_t x2 = a07;
    sp_t x3 = a10;

    sp_t t0 = sp_ntt_add_partial(x0, x2, p);
    sp_t t2 = sp_ntt_sub_partial(x0, x2, p);
    sp_t t1 = sp_ntt_add_partial(x1, x3, p);
    sp_t t3 = sp_ntt_sub_partial(x1, x3, p);

    t0 = sp_ntt_mul(t0, ntt_const[4], ntt_const[NC+4], p);
    t2 = sp_ntt_mul(t2, ntt_const[5], ntt_const[NC+5], p);
    t1 = sp_ntt_mul(t1, ntt_const[6], ntt_const[NC+6], p);
    t3 = sp_ntt_mul(t3, ntt_const[7], ntt_const[NC+7], p);

    x0 = sp_ntt_add(t0, t1, p);
    x1 = sp_ntt_sub(t0, t1, p);
    x2 = sp_ntt_add(t2, t3, p);
    x3 = sp_ntt_sub(t2, t3, p);

    a01 = x0;
    a04 = x2;
    a07 = x1;
    a10 = x3;
  }
  {
    sp_t x0 = a02;
    sp_t x1 = a05;
    sp_t x2 = a08;
    sp_t x3 = a11;

    sp_t t0 = sp_ntt_add_partial(x0, x2, p);
    sp_t t2 = sp_ntt_sub_partial(x0, x2, p);
    sp_t t1 = sp_ntt_add_partial(x1, x3, p);
    sp_t t3 = sp_ntt_sub_partial(x1, x3, p);

    t0 = sp_ntt_mul(t0, ntt_const[8], ntt_const[NC+8], p);
    t2 = sp_ntt_mul(t2, ntt_const[9], ntt_const[NC+9], p);
    t1 = sp_ntt_mul(t1, ntt_const[10], ntt_const[NC+10], p);
    t3 = sp_ntt_mul(t3, ntt_const[11], ntt_const[NC+11], p);

    x0 = sp_ntt_add(t0, t1, p);
    x1 = sp_ntt_sub(t0, t1, p);
    x2 = sp_ntt_add(t2, t3, p);
    x3 = sp_ntt_sub(t2, t3, p);

    a02 = x0;
    a05 = x2;
    a08 = x1;
    a11 = x3;
  }
  {
    sp_t x00, x01, x02;

    a01 = sp_ntt_add(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add(a01, a02, p);
    x02 = sp_ntt_sub(a01, a02, p);

    out[ 0 * ostride] = x00;
    out[ 4 * ostride] = x01;
    out[ 8 * ostride] = x02;
  }
  {
    sp_t x03, x04, x05;

    a04 = sp_ntt_add(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add(a04, a05, p);
    x05 = sp_ntt_sub(a04, a05, p);

    out[ 1 * ostride] = x04;
    out[ 5 * ostride] = x05;
    out[ 9 * ostride] = x03;
  }
  {
    sp_t x06, x07, x08;

    a07 = sp_ntt_add(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add(a07, a08, p);
    x08 = sp_ntt_sub(a07, a08, p);

    out[ 2 * ostride] = x08;
    out[ 6 * ostride] = x06;
    out[10 * ostride] = x07;
  }
  {
    sp_t x09, x10, x11;

    a10 = sp_ntt_add(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add(a10, a11, p);
    x11 = sp_ntt_sub(a10, a11, p);

    out[ 3 * ostride] = x09;
    out[ 7 * ostride] = x10;
    out[11 * ostride] = x11;
  }
}

static void
ntt12_twiddle_run_core(spv_t in, spv_size_t istride,
                spv_t out, spv_size_t ostride,
                spv_t w, sp_t p, spv_t ntt_const)
{
  sp_t a00, a01, a02, a03, 
       a04, a05, a06, a07,
       a08, a09, a10, a11;

  {
    sp_t x00 = in[ 0 * istride];
    sp_t x01 = in[ 4 * istride];
    sp_t x02 = in[ 8 * istride];

    a01 = sp_ntt_add(x01, x02, p);
    a02 = sp_ntt_sub(x01, x02, p);

    a00 = sp_ntt_add(x00, a01, p);
  }
  {
    sp_t x03 = in[ 3 * istride];
    sp_t x04 = in[ 7 * istride];
    sp_t x05 = in[11 * istride];

    a04 = sp_ntt_add(x04, x05, p);
    a05 = sp_ntt_sub(x04, x05, p);

    a03 = sp_ntt_add(x03, a04, p);
  }
  {
    sp_t x08 = in[ 2 * istride];
    sp_t x06 = in[ 6 * istride];
    sp_t x07 = in[10 * istride];

    a07 = sp_ntt_add(x07, x08, p);
    a08 = sp_ntt_sub(x07, x08, p);

    a06 = sp_ntt_add(x06, a07, p);
  }
  {
    sp_t x10 = in[ 1 * istride];
    sp_t x11 = in[ 5 * istride];
    sp_t x09 = in[ 9 * istride];

    a10 = sp_ntt_add(x10, x11, p);
    a11 = sp_ntt_sub(x10, x11, p);

    a09 = sp_ntt_add(x09, a10, p);
  }
  {
    sp_t x0 = a00;
    sp_t x1 = a03;
    sp_t x2 = a06;
    sp_t x3 = a09;

    sp_t t0 = sp_ntt_add(x0, x2, p);
    sp_t t2 = sp_ntt_sub(x0, x2, p);
    sp_t t1 = sp_ntt_add(x1, x3, p);
    sp_t t3 = sp_ntt_sub_partial(x1, x3, p);

    t3 = sp_ntt_mul(t3, ntt_const[3], ntt_const[NC+3], p);

    x0 = sp_ntt_add(t0, t1, p);
    x1 = sp_ntt_sub(t0, t1, p);
    x2 = sp_ntt_add(t2, t3, p);
    x3 = sp_ntt_sub(t2, t3, p);

    a00 = x0;
    a03 = x2;
    a06 = x1;
    a09 = x3;
  }
  {
    sp_t x0 = a01;
    sp_t x1 = a04;
    sp_t x2 = a07;
    sp_t x3 = a10;

    sp_t t0 = sp_ntt_add_partial(x0, x2, p);
    sp_t t2 = sp_ntt_sub_partial(x0, x2, p);
    sp_t t1 = sp_ntt_add_partial(x1, x3, p);
    sp_t t3 = sp_ntt_sub_partial(x1, x3, p);

    t0 = sp_ntt_mul(t0, ntt_const[4], ntt_const[NC+4], p);
    t2 = sp_ntt_mul(t2, ntt_const[5], ntt_const[NC+5], p);
    t1 = sp_ntt_mul(t1, ntt_const[6], ntt_const[NC+6], p);
    t3 = sp_ntt_mul(t3, ntt_const[7], ntt_const[NC+7], p);

    x0 = sp_ntt_add(t0, t1, p);
    x1 = sp_ntt_sub(t0, t1, p);
    x2 = sp_ntt_add(t2, t3, p);
    x3 = sp_ntt_sub(t2, t3, p);

    a01 = x0;
    a04 = x2;
    a07 = x1;
    a10 = x3;
  }
  {
    sp_t x0 = a02;
    sp_t x1 = a05;
    sp_t x2 = a08;
    sp_t x3 = a11;

    sp_t t0 = sp_ntt_add_partial(x0, x2, p);
    sp_t t2 = sp_ntt_sub_partial(x0, x2, p);
    sp_t t1 = sp_ntt_add_partial(x1, x3, p);
    sp_t t3 = sp_ntt_sub_partial(x1, x3, p);

    t0 = sp_ntt_mul(t0, ntt_const[8], ntt_const[NC+8], p);
    t2 = sp_ntt_mul(t2, ntt_const[9], ntt_const[NC+9], p);
    t1 = sp_ntt_mul(t1, ntt_const[10], ntt_const[NC+10], p);
    t3 = sp_ntt_mul(t3, ntt_const[11], ntt_const[NC+11], p);

    x0 = sp_ntt_add(t0, t1, p);
    x1 = sp_ntt_sub(t0, t1, p);
    x2 = sp_ntt_add(t2, t3, p);
    x3 = sp_ntt_sub(t2, t3, p);

    a02 = x0;
    a05 = x2;
    a08 = x1;
    a11 = x3;
  }
  {
    sp_t x00, x01, x02;

    a01 = sp_ntt_add(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add_partial(a01, a02, p);
    x02 = sp_ntt_sub_partial(a01, a02, p);

    x01 = sp_ntt_mul(x01, w[6], w[7], p);
    x02 = sp_ntt_mul(x02, w[14], w[15], p);

    out[ 0 * ostride] = x00;
    out[ 4 * ostride] = x01;
    out[ 8 * ostride] = x02;
  }
  {
    sp_t x03, x04, x05;

    a04 = sp_ntt_add(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add_partial(a04, a05, p);
    x05 = sp_ntt_sub_partial(a04, a05, p);

    x04 = sp_ntt_mul(x04, w[0], w[1], p);
    x05 = sp_ntt_mul(x05, w[8], w[9], p);
    x03 = sp_ntt_mul(x03, w[16], w[17], p);

    out[ 1 * ostride] = x04;
    out[ 5 * ostride] = x05;
    out[ 9 * ostride] = x03;
  }
  {
    sp_t x06, x07, x08;

    a07 = sp_ntt_add(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add_partial(a07, a08, p);
    x08 = sp_ntt_sub_partial(a07, a08, p);

    x08 = sp_ntt_mul(x08, w[2], w[3], p);
    x06 = sp_ntt_mul(x06, w[10], w[11], p);
    x07 = sp_ntt_mul(x07, w[18], w[19], p);

    out[ 2 * ostride] = x08;
    out[ 6 * ostride] = x06;
    out[10 * ostride] = x07;
  }
  {
    sp_t x09, x10, x11;

    a10 = sp_ntt_add(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add_partial(a10, a11, p);
    x11 = sp_ntt_sub_partial(a10, a11, p);

    x09 = sp_ntt_mul(x09, w[4], w[5], p);
    x10 = sp_ntt_mul(x10, w[12], w[13], p);
    x11 = sp_ntt_mul(x11, w[20], w[21], p);

    out[ 3 * ostride] = x09;
    out[ 7 * ostride] = x10;
    out[11 * ostride] = x11;
  }
}

static void
ntt12_pfa_run_core(spv_t x, spv_size_t start,
          spv_size_t inc, spv_size_t n,
          sp_t p, spv_t ntt_const)
{
  spv_size_t j00, j01, j02, j03, 
             j04, j05, j06, j07, 
             j08, j09, j10, j11;

  sp_t a00, a01, a02, a03, 
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
    sp_t x00 = x[j00];
    sp_t x01 = x[j04];
    sp_t x02 = x[j08];

    a01 = sp_ntt_add(x01, x02, p);
    a02 = sp_ntt_sub(x01, x02, p);

    a00 = sp_ntt_add(x00, a01, p);
  }
  {
    sp_t x03 = x[j03];
    sp_t x04 = x[j07];
    sp_t x05 = x[j11];

    a04 = sp_ntt_add(x04, x05, p);
    a05 = sp_ntt_sub(x04, x05, p);

    a03 = sp_ntt_add(x03, a04, p);
  }
  {
    sp_t x08 = x[j02];
    sp_t x06 = x[j06];
    sp_t x07 = x[j10];

    a07 = sp_ntt_add(x07, x08, p);
    a08 = sp_ntt_sub(x07, x08, p);

    a06 = sp_ntt_add(x06, a07, p);
  }
  {
    sp_t x10 = x[j01];
    sp_t x11 = x[j05];
    sp_t x09 = x[j09];

    a10 = sp_ntt_add(x10, x11, p);
    a11 = sp_ntt_sub(x10, x11, p);

    a09 = sp_ntt_add(x09, a10, p);
  }
  {
    sp_t x0 = a00;
    sp_t x1 = a03;
    sp_t x2 = a06;
    sp_t x3 = a09;

    sp_t t0 = sp_ntt_add(x0, x2, p);
    sp_t t2 = sp_ntt_sub(x0, x2, p);
    sp_t t1 = sp_ntt_add(x1, x3, p);
    sp_t t3 = sp_ntt_sub_partial(x1, x3, p);

    t3 = sp_ntt_mul(t3, ntt_const[3], ntt_const[NC+3], p);

    x0 = sp_ntt_add(t0, t1, p);
    x1 = sp_ntt_sub(t0, t1, p);
    x2 = sp_ntt_add(t2, t3, p);
    x3 = sp_ntt_sub(t2, t3, p);

    a00 = x0;
    a03 = x2;
    a06 = x1;
    a09 = x3;
  }
  {
    sp_t x0 = a01;
    sp_t x1 = a04;
    sp_t x2 = a07;
    sp_t x3 = a10;

    sp_t t0 = sp_ntt_add_partial(x0, x2, p);
    sp_t t2 = sp_ntt_sub_partial(x0, x2, p);
    sp_t t1 = sp_ntt_add_partial(x1, x3, p);
    sp_t t3 = sp_ntt_sub_partial(x1, x3, p);

    t0 = sp_ntt_mul(t0, ntt_const[4], ntt_const[NC+4], p);
    t2 = sp_ntt_mul(t2, ntt_const[5], ntt_const[NC+5], p);
    t1 = sp_ntt_mul(t1, ntt_const[6], ntt_const[NC+6], p);
    t3 = sp_ntt_mul(t3, ntt_const[7], ntt_const[NC+7], p);

    x0 = sp_ntt_add(t0, t1, p);
    x1 = sp_ntt_sub(t0, t1, p);
    x2 = sp_ntt_add(t2, t3, p);
    x3 = sp_ntt_sub(t2, t3, p);

    a01 = x0;
    a04 = x2;
    a07 = x1;
    a10 = x3;
  }
  {
    sp_t x0 = a02;
    sp_t x1 = a05;
    sp_t x2 = a08;
    sp_t x3 = a11;

    sp_t t0 = sp_ntt_add_partial(x0, x2, p);
    sp_t t2 = sp_ntt_sub_partial(x0, x2, p);
    sp_t t1 = sp_ntt_add_partial(x1, x3, p);
    sp_t t3 = sp_ntt_sub_partial(x1, x3, p);

    t0 = sp_ntt_mul(t0, ntt_const[8], ntt_const[NC+8], p);
    t2 = sp_ntt_mul(t2, ntt_const[9], ntt_const[NC+9], p);
    t1 = sp_ntt_mul(t1, ntt_const[10], ntt_const[NC+10], p);
    t3 = sp_ntt_mul(t3, ntt_const[11], ntt_const[NC+11], p);

    x0 = sp_ntt_add(t0, t1, p);
    x1 = sp_ntt_sub(t0, t1, p);
    x2 = sp_ntt_add(t2, t3, p);
    x3 = sp_ntt_sub(t2, t3, p);

    a02 = x0;
    a05 = x2;
    a08 = x1;
    a11 = x3;
  }
  {
    sp_t x00, x01, x02;

    a01 = sp_ntt_add(a00, a01, p);

    x00 = a00;
    x01 = sp_ntt_add(a01, a02, p);
    x02 = sp_ntt_sub(a01, a02, p);

    x[j00] = x00;
    x[j04] = x01;
    x[j08] = x02;
  }
  {
    sp_t x03, x04, x05;

    a04 = sp_ntt_add(a03, a04, p);

    x03 = a03;
    x04 = sp_ntt_add(a04, a05, p);
    x05 = sp_ntt_sub(a04, a05, p);

    x[j01] = x04;
    x[j05] = x05;
    x[j09] = x03;
  }
  {
    sp_t x06, x07, x08;

    a07 = sp_ntt_add(a06, a07, p);

    x06 = a06;
    x07 = sp_ntt_add(a07, a08, p);
    x08 = sp_ntt_sub(a07, a08, p);

    x[j02] = x08;
    x[j06] = x06;
    x[j10] = x07;
  }
  {
    sp_t x09, x10, x11;

    a10 = sp_ntt_add(a09, a10, p);

    x09 = a09;
    x10 = sp_ntt_add(a10, a11, p);
    x11 = sp_ntt_sub(a10, a11, p);

    x[j03] = x09;
    x[j07] = x10;
    x[j11] = x11;
  }
}

DECLARE_CORE_ROUTINES(12)
