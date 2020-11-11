/******************************************************************************/

/* prime64 : deterministic M-R primality test for a 64-bit value. */

/* Copyright (c) Brett Hale 2020.
 * distributed under BSD 2-clause license terms. */

/* note: requires gcc / clang '__int128' extended type. */

/******************************************************************************/

#include <inttypes.h>
#include <stdio.h>


/* return (1) if the nul-terminated C string forms a valid
 * 64-bit unsigned integer value in C locale decimal format,
 * and store the value in (u); return (0) otherwise: */

static int u64_arg (uint64_t *u, const char *s)
{
    int ret;

    if ((ret = *s) != 0)
    {
        uint64_t x = 0, d;

        if (ret == '0') /* "0" or not a decimal format: */
            return (s[1] ? (0) : (*u = x) == 0);

        for (; (d = (uint64_t) (*s++)) != 0; x += d)
        {
            if ((d -= ('0')) > (9) ||
                (x > (UINT64_MAX / (10)))) return (0);
            if ((x *= (10)) > (UINT64_MAX - d))
                return (0);
        }

        *u = x; /* a valid 64-bit unsigned integer value. */
    }

    return ret;
}

/******************************************************************************/

static int sprp (uint64_t n, uint64_t a)
{
    uint64_t m = n - 1, r, y;
    unsigned int s = 1, j;

    /* assert(n > 2 && (n & 0x1) != 0); */
    /* note: modified M-R test for successive bases. */

    while ((m & (UINT64_C(1) << s)) == 0) s++;
    r = m >> s; /* r, s s.t. 2^s * r = n - 1, r in odd. */

    if ((a %= n) == 0) /* else (0 < a < n) */
        return (1);

    {
        __extension__ unsigned __int128 u = 1, w = a;

        while (r != 0)
        {
            if ((r & 0x1) != 0)
                u = (u * w) % n; /* (mul-rdx) */

            if ((r >>= 1) != 0)
                w = (w * w) % n; /* (sqr-rdx) */
        }

        if ((y = (uint64_t) u) == 1)
            return (1);
    }

    for (j = 1; j < s && y != m; j++)
    {
        __extension__ unsigned __int128 u = y;
        u = (u * u) % n; /* (sqr-rdx) */

        if ((y = (uint64_t) u) <= 1) /* (n) is composite: */
            return (0);
    }

    return (y == m);
}

/******************************************************************************/

static int is_prime (uint64_t n)
{
    const uint32_t sprp32_base[] = /* (Jaeschke) */ {
        2, 7, 61, 0};

    const uint32_t sprp64_base[] = /* (Sinclair) */ {
        2, 325, 9375, 28178, 450775, 9780504, 1795265022, 0};

    const uint32_t *sprp_base;

    /* assert(n > 1); */

    if ((n & 0x1) == 0) /* even: */
        return (n == 2);

    sprp_base = (n <= UINT32_MAX) ? sprp32_base : sprp64_base;

    for (; *sprp_base != 0; sprp_base++)
        if (!sprp(n, *sprp_base)) return (0);

    return (1);
}


int main (int argc, char **argv)
{
    uint64_t n = 0;

    if (argc < 2 || !u64_arg(& n, argv[1]) || (n < 2))
    {
        fprintf(stderr, "usage: prime64 < u64 = 2 .. 2^64 - 1 >\n");
        return (1);
    }

    fprintf(stdout, "%"PRIu64" : %s\n", n,
            is_prime(n) ? "prime" : "composite");

    return (0);
}

/******************************************************************************/
