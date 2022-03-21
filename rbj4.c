/******************************************************************************/

/* RBJ.4 : exact p(k, 1) values for 2 <= k <= 24 (Monier's result) : */

/* Copyright (c) 2020 Brett Hale.
 * distributed under BSD-2-Clause license terms. see: mrtab.c */

/******************************************************************************/

#include <inttypes.h>
#include <stdio.h>

#include <float.h>
#include <math.h>

#include "spk12.h" /* small prime factorization. */

#if defined (QUADMATH)
#pragma GCC diagnostic ignored "-Wpedantic" /* (Q-suffix) */
#include <quadmath.h>
#endif


/* return gcd(u, v), where: (0 <= u, v < 2^32) : */

static uint32_t ugcd (uint32_t u, uint32_t v)
{
    /* additional iteration if (u < v) : */
    for (uint32_t t = 0; (t = v) != 0; u = t)
        v = u % v;

    return u;
}

/******************************************************************************/

static uint32_t sprp_bases (uint32_t n)
{
    uint32_t pbuf[(24)], pn, sn, un, vn, wn, p, i;

    /* assert(n > 1 && n < (UINT32_C(1) << (24))); */
    /* assert((n & 0x1) != 0); */

    /* a return value of (0) indicates that (n) is a prime, therefore
     * S(n) / (n - 1) should not be added to the composite running sum.
     * for a prime (p) : S(p) / (p - 1) = (1) */

    if ((pn = (uint32_t) sp_factor(pbuf, n)) == 1)
        return (0);

    /* prime vs. distinct prime factorization: */

    for (wn = 1, p = pbuf[0], i = 1; i < pn; i++)
    {
        uint32_t pi = pbuf[i];
        if (pi != p) pbuf[wn++] = (p = pi);
    }

    /* Monier's formula for S(n) : */

    for (vn = (24), i = 0; vn > 1 && i < wn; i++)
    {
        uint32_t pi = pbuf[i], vi = 0;
        do pi >>= 1, vi++; while ((pi & 0x1) == 0);
        if (vi < vn) vn = vi;
    }

    sn = UINT32_C(1);
    sn = 1 + ((sn << (wn * vn)) - 1) / ((sn << (wn)) - 1);

    /* find u(n), the largest odd factor of (n - 1). find the product
     * of gcd(p - 1, u(n)) over the unique prime factors (p) of (n). */

    for (un = n >> 1; (un & 0x1) == 0; un >>= 1);
    for (i = 0; i < wn; i++) sn *= ugcd(un, pbuf[i] - 1);

    return sn;
}

/******************************************************************************/

int main (void)
{
    unsigned int kmax = (24), k; /* {0 .. kmax} table: */
    double pk[(24) + 1];

    /* bias summation terms such that: fp{p(k, 1)} >= p(k, 1) */

    /* Burthe [2] mentions 'several hours' on a SPARC I. at this time,
     * a modest 2.66 GHz Core 2 Duo takes less than (6) seconds. */

    pk[0] = 1.0, pk[1] = 1.0; /* all fail. */
    pk[2] = 0.0, pk[3] = 0.0; /* all pass. */

    for (k = 4; k <= kmax; k++)
    {
        uint32_t nmax = (UINT32_C(1) << k), n, sn;
        double num, den, en, ed, an, x, t;

        num = den = 0.0, en = ed = 0.0; /* 2Sum series: */

        for (n = (nmax >> 1) + 1; n < nmax; n += 2)
        {
            if ((sn = sprp_bases(n)) != 0)
            {
                an = (double) sn / (double) (n - 1);
                an = nextafter(an, DBL_MAX);

                x = num + an; t = x - num;
                en += num - (x - t) + (an - t); num = x;
            }
            else
                an = 1.0;

            x = den + an; t = x - den;
            ed += den - (x - t) + (an - t); den = x;
        }

        pk[k] = nextafter((num + en) / (den + ed), DBL_MAX);
        fprintf(stdout, "%2u : %.16e\n", k, pk[k]);
    }

#if defined (QUADMATH)

    /* evaluate with quad precision to show that:
     * 0 <= (fp{p(k, 1)} - p(k, 1)) / p(k, 1) < (2.0) * (EPS) */

    __float128 pkq[(24) + 1];
    char buf[128];

    pkq[0] = 1.0, pkq[1] = 1.0; /* all fail. */
    pkq[2] = 0.0, pkq[3] = 0.0; /* all pass. */

    for (k = 4; k <= kmax; k++)
    {
        uint32_t nmax = (UINT32_C(1) << k), n, sn;
        __float128 num, den, en, ed, an, x, t;

        num = den = 0.0, en = ed = 0.0; /* 2Sum series: */

        for (n = (nmax >> 1) + 1; n < nmax; n += 2)
        {
            if ((sn = sprp_bases(n)) != 0)
            {
                an = (__float128) sn / (__float128) (n - 1);

                x = num + an; t = x - num;
                en += num - (x - t) + (an - t); num = x;
            }
            else
                an = 1.0;

            x = den + an; t = x - den;
            ed += den - (x - t) + (an - t); den = x;
        }

        pkq[k] = (num + en) / (den + ed);
        quadmath_snprintf(buf, sizeof(buf), "%.35Qe", pkq[k]);
        fprintf(stdout, "%2u : %s\n", k, buf);
    }

    for (k = 4; k <= kmax; k++)
    {
        const char *const fmt = "k = %2u : rel = %s, "
            "rel / DBL_EPS = %s\n";

        __float128 rerr;
        char rbuf[128], ebuf[128];

        rerr = (pk[k] - pkq[k]) / pkq[k];
        quadmath_snprintf(rbuf, sizeof(rbuf), "%.16Qe", rerr);

        rerr = rerr / DBL_EPSILON;
        quadmath_snprintf(ebuf, sizeof(ebuf), "%2.2Qf", rerr);

        fprintf(stdout, fmt, k, rbuf, ebuf);
    }

#endif

    fprintf(stdout, "\n    %.16e", pk[0]);

    for (k = 1; k <= kmax; k++)
        fprintf(stdout, (k % 3) ? ", %.16e" : ",\n    %.16e", pk[k]);

    fprintf(stdout, "\n\n");

    /* since p(k, 1) < 1/5 for 2 <= k <= 24, the Monier-Rabin theorem
     * yields: p(k, t) <= 4^(1-t) * p(k, 1) / (1 - p(k, 1)) < (4^-t). */

    return (0);
}

/******************************************************************************/
