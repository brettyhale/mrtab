/******************************************************************************/

/* RBJ.2.L2 : c(s) sequence: */

/* Copyright (c) 2020 Brett Hale.
 * distributed under BSD-2-Clause license terms. see: mrtab.c */

/******************************************************************************/

#include <inttypes.h>
#include <stdio.h>

#include <float.h>
#include <math.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

#if defined (QUADMATH)
#pragma GCC diagnostic ignored "-Wpedantic" /* (Q-suffix) */
#include <quadmath.h>
#endif


int main (void)
{
    unsigned int smax = (30), s; /* {0 .. smax} table: */
    double c[(30) + 1], r, e;

    /* bias summation terms such that: fp{c(s)} >= c(s) */

    /* TODO: not satisfied with this summation. is there a
     * formula to prevent (s) scaling the relative error? */

    r = M_PI * M_PI / 6.0, e = 0.0;
    c[0] = nextafter(r, DBL_MAX);

    for (s = 1; s <= smax; s++)
    {
        double y, t; /* Kahan summation: */

        y = - 1.0 / (double) (s * s);
        y = nextafter(y, DBL_MAX) - e;
        t = r + y; e = (t - r) - y; r = t;

        c[s] = (s + 1) * r;

        /* recurrence relation: */
        /* c[s] = c[s - 1] + r - (1.0 / s); */
    }

#if defined (QUADMATH)

    /* evaluate with quad precision to show that:
     * 0 <= (fp{c(s)} - c(s)) / c(s) < (s + 1) * (EPS) */

    __float128 cq[(30) + 1], rq, eq;

    rq = M_PIq * M_PIq / 6.0, eq = 0.0;
    cq[0] = nextafterq(rq, FLT128_MAX);

    for (s = 1; s <= smax; s++)
    {
        __float128 y, t; /* Kahan summation: */

        y = - 1.0 / (__float128) (s * s);
        y = nextafterq(y, FLT128_MAX) - eq;
        t = rq + y; eq = (t - rq) - y; rq = t;

        cq[s] = (s + 1) * rq;
    }

    for (s = 0; s <= smax; s++)
    {
        const char *const fmt = "s = %2u : rel = %s, "
            "rel / DBL_EPS = %s\n";

        __float128 rerr;
        char rbuf[128], ebuf[128];

        rerr = (c[s] - cq[s]) / cq[s];
        quadmath_snprintf(rbuf, sizeof(rbuf), "%.16Qe", rerr);

        rerr = rerr / DBL_EPSILON;
        quadmath_snprintf(ebuf, sizeof(ebuf), "%2.2Qf", rerr);

        fprintf(stdout, fmt, s, rbuf, ebuf);
    }

#endif

    fprintf(stdout, "\n    %.16e", c[0]);

    for (s = 1; s <= smax; s++)
        fprintf(stdout, (s % 3) ? ", %.16e" : ",\n    %.16e", c[s]);

    fprintf(stdout, "\n\n");

    return (0);
}

/******************************************************************************/
