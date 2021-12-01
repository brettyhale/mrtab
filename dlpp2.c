/******************************************************************************/

/* DLP.4.proposition.2 : test k-bit prime probability bounds: */

/* Copyright (c) 2020 Brett Hale.
 * distributed under BSD-2-Clause license terms. see: mrtab.c */

/******************************************************************************/

#include <inttypes.h>
#include <stdio.h>


int main (void)
{
    for (unsigned int k = 4; k <= (20); k++)
    {
        uint32_t nmax = (UINT32_C(1) << k), n, p;

        for (p = 0, n = (nmax >> 1) + 1; n < nmax; n += 2)
        {
            uint32_t d, q, c = 0;

            for (d = 3; !c && (q = n / d) >= d; d += 2)
                c = (q * d == n);

            p += (c == 0); /* (n) is prime. */
        }

        double lhs = (double) p * k, rhs = 0.71867 * nmax;
        fprintf(stdout, "%2u : %s\n", k, ((lhs > rhs) ? "T" : "F"));
    }

    return (0);
}

/******************************************************************************/
