/******************************************************************************/

/* mrtab : Miller-Rabin test iterations for a random, k-bit probable prime
 * search, satisfying an upper bound for the error probability: p(k, t)
 *
 * Copyright (c) 2020 Brett Hale.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */


/* [1] Damgaard, Landrock, Pomerance, "Average Case Error Estimates for
 * the Strong Probable Prime Test". Mathematics of Computation, Vol. 61,
 * Jul. 1993, pp. 177-194. aka: [DLP] (referenced in HAC.4.4.1). */

/* M-R iterations for (k > 16) are computed, so that trial division can
 * be performed with 16-bit arithmetic. the probablility that a random,
 * k-bit odd is prime, according to DLP.4.proposition.2, satisfies:

 * pi(2^k) - pi(2^(k - 1)) > (0.71867) * (2^k) / k, for all k >= 21.

 * by counting the primes < (2^20), it is clear that this lower bound is
 * actually valid for all k >= 8. */


/* [2] R. Burthe, Jr., "Further Investigations with the Strong Probable
 * Prime Test". Mathematics of Computation, Vol. 65, Jan. 1996, pp. 373-
 * 381. aka: [RBJ] */

/* an important result in [2] is that it validates the conjecture in [1]
 * that: p(k, t) <= (4^-t), for all k >= 2, t >= 1. */


/* [HAC] Menezes, Oorschot, Vanstone, "Handbook of Applied Cryptography".
 * CRC Press, ISBN 0849385237. also available online (conditionally) at:
 * [http://www.cacr.math.uwaterloo.ca/hac/] */

/* chapter 4, in particular section 4.4.1, introduces the random probable
 * prime search using the Miller-Rabin test. */

/******************************************************************************/

#include <string.h>
#include <stdio.h>

#include <float.h>
#include <math.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif


/* return (1) if the nul-terminated C string forms a valid
 * 32-bit unsigned integer value in C locale decimal format,
 * and store the value in (u); return (0) otherwise: */

static int u32_arg (unsigned long *u, const char *s)
{
    int ret;

    if ((ret = *s) != 0)
    {
        unsigned long x = 0, d;

        if (ret == '0') /* "0" or not a decimal format: */
            return (s[1] ? (0) : (*u = x) == 0);

        for (; (d = (unsigned long) (*s++)) != 0; x += d)
        {
            if ((d -= ('0')) > (9) ||
                (x > (0xffffffffUL / (10)))) return (0);
            if ((x *= (10)) > (0xffffffffUL - d))
                return (0);
        }

        *u = x; /* a valid 32-bit unsigned integer value. */
    }

    return ret;
}

/******************************************************************************/

/* the signature for the p(k, t) evaluation function: */

/* note: an implementation may assert that: (k > 16 && t >= 1), although
 * it should attempt to handle all (k > 1). the values of (k <= 1) yield:
 * p(k, t) = (1) (no primes), while (t = 0) is meaningless. */

typedef double (*p_kt_fn)(unsigned int, unsigned int);

static int dlp_tab (p_kt_fn p_kt)
{
    unsigned int k, t;

    fprintf(stdout, "lower bounds for -lb(p(k, t))\n\n");
    /* DLP.table.1 actually lists: floor(-lb(p(k, t))) */

    fprintf(stdout, "k\\t |");
    for (t = 1; t <= 10; t++)
        fprintf(stdout, " %3u", t);
    fprintf(stdout, "\n-----");
    for (t = 1; t <= 10; t++)
        fprintf(stdout, "----");
    fprintf(stdout, "\n");

    for (k = 100; k <= 600; k += 50)
    {
        fprintf(stdout, "%3u |", k);
        for (t = 1; t <= 10; t++)
            fprintf(stdout, " %3u", /* floor(-lb(p(k, t))) : */
                    (unsigned int) (- log2((*p_kt)(k, t))));
        fprintf(stdout, "\n");
    }

    return (0);
}


/* RBJ.4 : exact p(k, 1) values for 2 <= k <= 24 (Monier's result) : */

/* since p(k, 1) < 1/5 for 2 <= k <= 24, the Monier-Rabin theorem yields:
 * p(k, t) <= 4^(1-t) * p(k, 1) / (1 - p(k, 1)) < (4^-t). */

static const double p_k1_lut[25] = /* exact p(0, 1) .. p(24, 1) : */
{
    1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
    0.0000000000000000e+00, 1.6417910447761200e-01, 6.4299424184261059e-02,
    6.5348064836078495e-02, 5.6654752251003034e-02, 3.8003778178391873e-02,
    3.0837119635400381e-02, 2.0525079764265652e-02, 1.7393574680619316e-02,
    1.0710359182783314e-02, 7.9490871698650184e-03, 5.9337043808932611e-03,
    3.9442643069209568e-03, 2.6255166117476652e-03, 1.9286518790611249e-03,
    1.2577894174913744e-03, 9.0457147250914852e-04, 6.0885312016630043e-04,
    4.0170629568174411e-04, 2.7576379216948154e-04, 1.8760654682551843e-04,
    1.2612847365349537e-04
};

/******************************************************************************/

#if (0) /* RBJ.3,4 estimate: */

static const double rbj_lut[31] = /* RBJ.2.L2 : c(s) sequence: */
{
    1.6449340668482266e+00, 1.2898681336964530e+00, 1.1848022005446794e+00,
    1.1352918229484619e+00, 1.1066147786855773e+00, 1.0879377344226930e+00,
    1.0748162457153638e+00, 1.0650961175522526e+00, 1.0576081322462842e+00,
    1.0516633568168590e+00, 1.0468296924985450e+00, 1.0428224744612227e+00,
    1.0394465695552135e+00, 1.0365637612961471e+00, 1.0340734177152595e+00,
    1.0319005344518326e+00, 1.0299880678550721e+00, 1.0282918642340901e+00,
    1.0267772147162311e+00, 1.0254164587040657e+00, 1.0241872816392692e+00,
    1.0230714832592795e+00, 1.0220540713413131e+00, 1.0211225848400847e+00,
    1.0202665814306437e+00, 1.0194772446878697e+00, 1.0187470795427285e+00,
    1.0180696737096060e+00, 1.0174395089951531e+00, 1.0168518107322035e+00,
    1.0163024266454992e+00
};


/* the inner-most summation of 'N1' : note that c(s) evaluation has been
 * moved inside the summation over (m). this requires very little extra
 * effort relative to the exp2 function overhead, as (s) is clamped, and
 * the c(s) values have been precomputed. */

static double rbj_ktm (double rk, double rt, double rm)
{
    double rj, rs, r0 = 0.0;
    unsigned int jh, ji;

    rs = 30.0; /* RBJ.3 default. */

    jh = (unsigned int) ceil(rm);
    for (rj = 2.0, ji = 2; ji <= jh; rj += 1.0, ji++)
    {
        double jd, jn, js;

        jd = exp2((rk - 1.0) / rj) - 1.0;
        jn = exp2(rj - rm - 2.0);

        if ((js = jd * jn) < rs)
            rs = js; /* maximum (s) */

        jn = ceil(1.0 / (2.0 * jn)) - 1.0;
        r0 += jn / jd;
    }

    r0 *= exp2(- rm * rt);

    return (r0 * rbj_lut[(unsigned int) rs]); /* c(s) */
}


/* rbj_kt: this is an implementation of the estimate in RBJ.3, using the
 * default value for (q). the results match those in RBJ.table.4 - except
 * for where the combined (italicized) values have been used. */

static double rbj_kt (unsigned int k, unsigned int t)
{
    double rk = k, rt = t, rp, rq;
    unsigned int q = 4, mh, mi;

    /* q = q ? q : (4); */ /* RBJ.3 default. */

    /* assert(k > 1 && t >= 1); */

    rp = exp2(- 2.0 * rt); /* (4^-t) [RBJ] */

    if (k < 25) /* Monier-Rabin: */
    {
        double p_k1 = p_k1_lut[k];

        if (t > 1)
            rp *= 4.0 * p_k1 / (1.0 - p_k1);
        else
            rp = p_k1;

        if (k < 10) return rp; /* no RBJ.3 result. */
    }

    rq = 1.0 / q; /* (fractional step) */

    /* integral 'M' candidates: */

    mh = (unsigned int) (2.0 * sqrt(rk - 1.0) - 3.0);
    for (mi = 3; mi <= mh; mi++)
    {
        /* (q(M - 2)(M + 1)/2) summation terms: */

        double rm, p1, n1 = 0.0;
        unsigned int qi;

        /* fractional summation: */
        for (qi = q * 2 + 1; qi <= q * mi; qi++)
        {
            /* ensure integral (m) values are exact: */
            rm = ((qi % q) == 0) ? (qi / q) : (rq * qi);
            n1 += rbj_ktm(rk, rt, rm);
        }

        n1 *= 0.5 * (exp2(rt * rq) - 1.0);
        n1 += exp2(- (rt * mi + 2.0));

        p1 = 0.71867 / rk;
        if ((n1 /= (n1 + p1)) < rp) /* new 'M' candidate: */
            rp = n1;
    }

    return rp; /* p(k, t) */
}


/* finding asymptotically optimal values for (q) and (s) is impractical.
 * Burthe suggests the default values as reasonable choices, beyond which
 * the extra effort yields no significant improvement in the estimate. */

#endif /* (RBJ.3,4 estimate) */

/******************************************************************************/

/* dlp_kt: this is an implementation of the estimate in DLP.4, with a few
 * simple optimizations. the results match those in DLP.table.1. */


static double dlp_kt (unsigned int k, unsigned int t)
{
    const double c = (8.0 * (M_PI * M_PI - 6.0) / 3.0);

    double rk = k, rt = t, rp, mt;
    unsigned int mh, mi;

    /* assert(k > 1 && t >= 1); */

    rp = exp2(- 2.0 * rt); /* (4^-t) [RBJ] */

    if (k < 25) /* Monier-Rabin: */
    {
        double p_k1 = p_k1_lut[k];

        if (t > 1)
            rp *= 4.0 * p_k1 / (1.0 - p_k1);
        else
            rp = p_k1;

        if (k < 8) return rp; /* no DLP.4 result. */
    }

    mt = exp2(1.0 - rt);

    /* integral 'M' candidates: */

    mh = (unsigned int) (2.0 * sqrt(rk - 1.0) - 1.0);
    for (mi = 3; mi <= mh; mi++)
    {
        /* ((M - 2)(M + 1)/2) summation terms: */

        double r0 = 0.0;
        unsigned int j;

        for (j = 2; j <= mi; j++) /* {j .. M} */
        {
            unsigned int m = j + (j == 2);
            double ej = j, mj;

            ej += (rk - 1.0) / ej;
            mj = exp2((1.0 - rt) * m - ej);

            for (r0 += mj; m < mi; m++) /* {m .. M} \ {2} */
                r0 += (mj *= mt);
        }

        r0 *= c / (2.0 * mt);
        r0 += exp2(- (2.0 + rt * mi));

        if ((r0 *= rk / 0.71867) < rp) /* new 'M' candidate: */
            rp = r0;
    }

    return rp; /* p(k, t) */
}


/* note: both DLP.7 and RBJ.5 describe combined results from:

 * [3] S.H. Kim and C. Pomerance, "The Probability that a Random Probable
 * Prime is Composite," Mathematics of Computation, Vol. 53, Oct. 1989,
 * pp 721-742.

 * the complexity of this estimate is daunting - finding an optimal value
 * is probably not feasible. */

/******************************************************************************/

/* how to use the M-R table: the first entry is the maximum (k) value for
 * which (2) iterations of the M-R test are required to guarantee that:
 * p(k, t) <= (2^-s). an entry with a LUT index of (i) is the maximum (k)
 * value for which (i + 2) iterations are required. */

/* the (0) entry is an 'end-of-table' marker which delimits the maximum
 * number of iterations required for (k > 16). that is, if the (0) entry
 * has a LUT index of (j), then the maximum number of iterations required
 * is (tmax = j + 1). usage: "for (t = 1; k <= lut[t - 1]; t++);" */


static const char *usage =
    "usage: mrtab [s], where: s = 64 .. 256 (default: 128)\n"
    "M-R test iterations s.t. p(k, t) <= (2^-s), for k > 16.\n";

int main (int argc, char **argv)
{
    p_kt_fn p_kt = dlp_kt; /* default p(k, t) evaluation function. */

    unsigned int s = (128), kmax, tmax, k, t, ttab[256];
    double pmax;

    if (argc > 1) /* exponent option: */
    {
        unsigned long u;

        if (strcmp(argv[1], "-d") == 0)
            return dlp_tab(p_kt);

        if (!u32_arg(& u, argv[1]) || (u < 64) || (u > 256))
        {
            fprintf(stderr, "%s", usage);
            return (1);
        }

        s = (unsigned int) u;
    }

    pmax = exp2(- (double) s);
    fprintf(stdout, "k from t = 2 (k > 16) s.t. "
            "p(k, t) <= 2^-%u (%.2e) :\n", s, pmax);


    /* find kmax s.t. p(kmax, 1) <= (2^-s). kmax must be greater than,
     * or equal to, the minimum k value s.t. p(k, 1) <= (2^-s) : */

    for (kmax = (25); (*p_kt)(kmax, 1) > pmax; kmax <<= 1);

    /* find tmax s.t. p(k, t) <= (2^-s), for k > 16. the result from
     * [RBJ] yields: tmax = ceil(s/2), for k >= 2. (s) is an integral
     * value in this context: */

    tmax = (s + 1) >> 1;

    /* under the assertion that: p(k + 1, t) < p(k, t), on which the
     * implementation of a threshold table is predicated, this is the
     * upper bound for t s.t. p(k, t) <= (2^-s), for k > 16. */


    for (t = 2; t <= tmax; t++)
    {
        unsigned int k0 = (16) + 1, k1 = kmax;
        int found = 0;

        /* if: p(k + 1, t - 1) <= pmax < p(k, t - 1), assert that:
         * p(k, t) < p(k, t - 1). it does not follow, although it is
         * almost certainly the case, that: p(k, t) <= pmax. */

        while (!found && k0 <= k1)
        {
            k = k0 + (k1 - k0) / 2;

            if ((*p_kt)(k + 1, t - 1) > pmax) /* (k > k0) */
                k0 = k + 1;

            else if ((*p_kt)(k, t - 1) <= pmax) /* (k < k1) */
                k1 = k - 1;

            else /* sweet spot: */
            {
                while ((*p_kt)(k, t) > pmax) k++;

                if (k > kmax) /* pathological case (?) */
                {
                    unsigned int ti;

                    /* warning: found a local maxima in p(k, t).
                     * ensure that the table contains a sequence of
                     * non-increasing (k) values: */

                    /* should this be treated as an error? */
                    for (ti = t - 1; ti >= 2 && ttab[ti] < k; ti--)
                        ttab[ti] = k;
                }

                ttab[t] = kmax = k, found = 1;
            }
        }

        /* if no threshold (k) value was found, then tmax has been
         * determined to be unnecessarily high by p(k, t). clamp the
         * table to its current size: */

        if (!found) tmax = t - 1;
    }

#if (1)
#define MRTAB_FMT "0x%04x" /* threshold value LUT: */
#else
#define MRTAB_FMT "%6u"
#endif

    fprintf(stdout, "\n    "MRTAB_FMT, ttab[2]);

    ttab[++tmax] = 0; /* EOT entry. */

    for (t = 3; t <= tmax; t++)
    {
        if ((t - 2) % 8 != 0)
            fprintf(stdout, ", "MRTAB_FMT, ttab[t]);
        else
            fprintf(stdout, ",\n    "MRTAB_FMT, ttab[t]);
    }

    fprintf(stdout, "\n\n");

    /* the M-R implementation must handle candidates with (16) or fewer
     * significant bits explicitly, requiring up to (54) trial divisions;
     * i.e., the number of primes less than (2^8). */

    return (0);
}

/******************************************************************************/
