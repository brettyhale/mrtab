/******************************************************************************/

/* small prime LUT generator for trial division [HAC.4.4.1] : */

/* Copyright (c) 2020 Brett Hale.
 * distributed under BSD-2-Clause license terms. see: mrtab.c */

/******************************************************************************/

#include <inttypes.h>
#include <stdio.h>


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

int main (int argc, char **argv)
{
    uint32_t k = (12), nmax, n, p; /* prime table < (2^k) */
    double m = 1.0;

    const char *pfmt;

    if (argc > 1) /* (k) option: */
    {
        unsigned long u;

        if (!u32_arg(& u, argv[1]) || (u < 2) || (u > 16))
        {
            fprintf(stderr, "usage: sptab [k], where: "
                    "k = 2 .. 16 (default: 12)\n");
            return (1);
        }

        k = (uint32_t) u;
    }

    nmax = (UINT32_C(1) << k);
    pfmt = (nmax <= 0x100) ? "0x%02"PRIx32 : "0x%04"PRIx32;

    fprintf(stdout, "\n    ");
    fprintf(stdout, pfmt, UINT32_C(2));

    for (p = 1, n = 3; n < nmax; n += 2)
    {
        uint32_t d, q, c = 0;

        for (d = 3; !c && (q = n / d) >= d; d += 2)
            c = (q * d == n);

        if (!c) /* prime: */
        {
            m *= (double) (n - 1) / (double) n;
            fprintf(stdout, (p++ % 8) ? ", " : ",\n    ");
            fprintf(stdout, pfmt, n);
        }
    }

    fprintf(stdout, "\n\n");

    m = (1.0 - m) * 100.0;
    fprintf(stdout, "%"PRIu32" primes < 2^%"PRIu32" ", p, k);
    fprintf(stdout, "factor %2.2f%% of all odd integers.\n", m);

    return (0);
}

/******************************************************************************/
