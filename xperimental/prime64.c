/******************************************************************************/

/* prime64 : deterministic M-R primality test for a 64-bit value. */

/* requires '__int128' extended type. */

/* Copyright (c) 2020 Brett Hale.
 * distributed under BSD-2-Clause license terms. see: mrtab.c */

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

/* 54 primes < 2^8 factor 79.93% of all odd integers: */

static const uint8_t sp_lut[] =
{
    0x02, 0x03, 0x05, 0x07, 0x0b, 0x0d, 0x11, 0x13,
    0x17, 0x1d, 0x1f, 0x25, 0x29, 0x2b, 0x2f, 0x35,
    0x3b, 0x3d, 0x43, 0x47, 0x49, 0x4f, 0x53, 0x59,
    0x61, 0x65, 0x67, 0x6b, 0x6d, 0x71, 0x7f, 0x83,
    0x89, 0x8b, 0x95, 0x97, 0x9d, 0xa3, 0xa7, 0xad,
    0xb3, 0xb5, 0xbf, 0xc1, 0xc5, 0xc7, 0xd3, 0xdf,
    0xe3, 0xe5, 0xe9, 0xef, 0xf1, 0xfb, 0x00
};


static int sp_test (uint16_t n)
{
    uint16_t sp = sp_lut[0], q;

    /* assert(n > 1); */

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

    if (n < 65536) /* trial division for n < (2^16) : */
        return sp_test((uint16_t) n);

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
