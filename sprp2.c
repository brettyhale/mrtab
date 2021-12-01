/******************************************************************************/

/* demonstration of (2) as an effective witness to compositeness: */

/* Copyright (c) 2020 Brett Hale.
 * distributed under BSD-2-Clause license terms. see: mrtab.c */

/******************************************************************************/

#include <inttypes.h>
#include <stdio.h>

#include "spk12.h" /* small prime factorization. */


static int is_prime (uint32_t n)
{
    uint32_t sp = sp_lut[0], q;

    /* assert(n > 1 && n < (UINT32_C(1) << (24))); */

    for (unsigned int i = 1; (q = n / sp) >= sp; )
    {
        if (q * sp == n) /* composite: */
            return (0);

        else if ((sp = sp_lut[i++]) == 0) /* EOT entry: */
            break;
    }

    return (1);
}

/******************************************************************************/

static int sprp (uint32_t n, uint32_t a)
{
    uint32_t m = n - 1, r, y;
    unsigned int s = 1, j;

    /* assert(n > 3 && (n & 0x1) != 0); */
    /* assert(a > 1 && a < (n - 1)); */

    while ((m & (UINT32_C(1) << s)) == 0) s++;
    r = m >> s; /* r, s s.t. 2^s * r = n - 1, r in odd. */

    {
        uint64_t u = 1, w = a;

        while (r != 0)
        {
            if ((r & 0x1) != 0)
                u = (u * w) % n; /* (mul-rdx) */

            if ((r >>= 1) != 0)
                w = (w * w) % n; /* (sqr-rdx) */
        }

        if ((y = (uint32_t) u) == 1)
            return (1);
    }

    for (j = 1; j < s && y != m; j++)
    {
        uint64_t u = y;
        u = (u * u) % n; /* (sqr-rdx) */

        if ((y = (uint32_t) u) <= 1) /* (n) is composite: */
            return (0);
    }

    return (y == m);
}

/******************************************************************************/

int main (void)
{
    /* these (anecdotal) results show that a 2-SPRP test eliminates
     * the overwhelming majority of composite candidates, prior to
     * the independent M-R trials with randomized (a-SPRP) bases. */

    fprintf(stdout, "frequency of 2-SPRP strong liars:\n\n");

    for (unsigned int k = 4; k <= (24); k++)
    {
        uint32_t nmax = (UINT32_C(1) << k), n, c, s;

        for (c = 0, s = 0, n = (nmax >> 1) + 1; n < nmax; n += 2)
        {
            if (!is_prime(n))
            {
                c++; /* (composite) */
                s += (sprp(n, 2) != 0); /* 2-SPRP strong liar. */
            }
        }

        fprintf(stdout, "%2u : %2"PRIu32" / %7"PRIu32"\n", k, s, c);
    }

    return (0);
}

/******************************************************************************/
