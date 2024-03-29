/******************************************************************************/

/* prime factorization for (n) in: [2, 2^24 - 1] (with multiplicity).
 * small prime LUT generated using the sptab utility. */

/* Copyright (c) 2020 Brett Hale.
 * distributed under BSD-2-Clause license terms. see: mrtab.c */

/******************************************************************************/

#ifndef SP_K12_H_
#define SP_K12_H_

#include <stdint.h>

/******************************************************************************/

static const uint16_t sp_lut[] =
{
    0x0002, 0x0003, 0x0005, 0x0007, 0x000b, 0x000d, 0x0011, 0x0013,
    0x0017, 0x001d, 0x001f, 0x0025, 0x0029, 0x002b, 0x002f, 0x0035,
    0x003b, 0x003d, 0x0043, 0x0047, 0x0049, 0x004f, 0x0053, 0x0059,
    0x0061, 0x0065, 0x0067, 0x006b, 0x006d, 0x0071, 0x007f, 0x0083,
    0x0089, 0x008b, 0x0095, 0x0097, 0x009d, 0x00a3, 0x00a7, 0x00ad,
    0x00b3, 0x00b5, 0x00bf, 0x00c1, 0x00c5, 0x00c7, 0x00d3, 0x00df,
    0x00e3, 0x00e5, 0x00e9, 0x00ef, 0x00f1, 0x00fb, 0x0101, 0x0107,
    0x010d, 0x010f, 0x0115, 0x0119, 0x011b, 0x0125, 0x0133, 0x0137,
    0x0139, 0x013d, 0x014b, 0x0151, 0x015b, 0x015d, 0x0161, 0x0167,
    0x016f, 0x0175, 0x017b, 0x017f, 0x0185, 0x018d, 0x0191, 0x0199,
    0x01a3, 0x01a5, 0x01af, 0x01b1, 0x01b7, 0x01bb, 0x01c1, 0x01c9,
    0x01cd, 0x01cf, 0x01d3, 0x01df, 0x01e7, 0x01eb, 0x01f3, 0x01f7,
    0x01fd, 0x0209, 0x020b, 0x021d, 0x0223, 0x022d, 0x0233, 0x0239,
    0x023b, 0x0241, 0x024b, 0x0251, 0x0257, 0x0259, 0x025f, 0x0265,
    0x0269, 0x026b, 0x0277, 0x0281, 0x0283, 0x0287, 0x028d, 0x0293,
    0x0295, 0x02a1, 0x02a5, 0x02ab, 0x02b3, 0x02bd, 0x02c5, 0x02cf,
    0x02d7, 0x02dd, 0x02e3, 0x02e7, 0x02ef, 0x02f5, 0x02f9, 0x0301,
    0x0305, 0x0313, 0x031d, 0x0329, 0x032b, 0x0335, 0x0337, 0x033b,
    0x033d, 0x0347, 0x0355, 0x0359, 0x035b, 0x035f, 0x036d, 0x0371,
    0x0373, 0x0377, 0x038b, 0x038f, 0x0397, 0x03a1, 0x03a9, 0x03ad,
    0x03b3, 0x03b9, 0x03c7, 0x03cb, 0x03d1, 0x03d7, 0x03df, 0x03e5,
    0x03f1, 0x03f5, 0x03fb, 0x03fd, 0x0407, 0x0409, 0x040f, 0x0419,
    0x041b, 0x0425, 0x0427, 0x042d, 0x043f, 0x0443, 0x0445, 0x0449,
    0x044f, 0x0455, 0x045d, 0x0463, 0x0469, 0x047f, 0x0481, 0x048b,
    0x0493, 0x049d, 0x04a3, 0x04a9, 0x04b1, 0x04bd, 0x04c1, 0x04c7,
    0x04cd, 0x04cf, 0x04d5, 0x04e1, 0x04eb, 0x04fd, 0x04ff, 0x0503,
    0x0509, 0x050b, 0x0511, 0x0515, 0x0517, 0x051b, 0x0527, 0x0529,
    0x052f, 0x0551, 0x0557, 0x055d, 0x0565, 0x0577, 0x0581, 0x058f,
    0x0593, 0x0595, 0x0599, 0x059f, 0x05a7, 0x05ab, 0x05ad, 0x05b3,
    0x05bf, 0x05c9, 0x05cb, 0x05cf, 0x05d1, 0x05d5, 0x05db, 0x05e7,
    0x05f3, 0x05fb, 0x0607, 0x060d, 0x0611, 0x0617, 0x061f, 0x0623,
    0x062b, 0x062f, 0x063d, 0x0641, 0x0647, 0x0649, 0x064d, 0x0653,
    0x0655, 0x065b, 0x0665, 0x0679, 0x067f, 0x0683, 0x0685, 0x069d,
    0x06a1, 0x06a3, 0x06ad, 0x06b9, 0x06bb, 0x06c5, 0x06cd, 0x06d3,
    0x06d9, 0x06df, 0x06f1, 0x06f7, 0x06fb, 0x06fd, 0x0709, 0x0713,
    0x071f, 0x0727, 0x0737, 0x0745, 0x074b, 0x074f, 0x0751, 0x0755,
    0x0757, 0x0761, 0x076d, 0x0773, 0x0779, 0x078b, 0x078d, 0x079d,
    0x079f, 0x07b5, 0x07bb, 0x07c3, 0x07c9, 0x07cd, 0x07cf, 0x07d3,
    0x07db, 0x07e1, 0x07eb, 0x07ed, 0x07f7, 0x0805, 0x080f, 0x0815,
    0x0821, 0x0823, 0x0827, 0x0829, 0x0833, 0x083f, 0x0841, 0x0851,
    0x0853, 0x0859, 0x085d, 0x085f, 0x0869, 0x0871, 0x0883, 0x089b,
    0x089f, 0x08a5, 0x08ad, 0x08bd, 0x08bf, 0x08c3, 0x08cb, 0x08db,
    0x08dd, 0x08e1, 0x08e9, 0x08ef, 0x08f5, 0x08f9, 0x0905, 0x0907,
    0x091d, 0x0923, 0x0925, 0x092b, 0x092f, 0x0935, 0x0943, 0x0949,
    0x094d, 0x094f, 0x0955, 0x0959, 0x095f, 0x096b, 0x0971, 0x0977,
    0x0985, 0x0989, 0x098f, 0x099b, 0x09a3, 0x09a9, 0x09ad, 0x09c7,
    0x09d9, 0x09e3, 0x09eb, 0x09ef, 0x09f5, 0x09f7, 0x09fd, 0x0a13,
    0x0a1f, 0x0a21, 0x0a31, 0x0a39, 0x0a3d, 0x0a49, 0x0a57, 0x0a61,
    0x0a63, 0x0a67, 0x0a6f, 0x0a75, 0x0a7b, 0x0a7f, 0x0a81, 0x0a85,
    0x0a8b, 0x0a93, 0x0a97, 0x0a99, 0x0a9f, 0x0aa9, 0x0aab, 0x0ab5,
    0x0abd, 0x0ac1, 0x0acf, 0x0ad9, 0x0ae5, 0x0ae7, 0x0aed, 0x0af1,
    0x0af3, 0x0b03, 0x0b11, 0x0b15, 0x0b1b, 0x0b23, 0x0b29, 0x0b2d,
    0x0b3f, 0x0b47, 0x0b51, 0x0b57, 0x0b5d, 0x0b65, 0x0b6f, 0x0b7b,
    0x0b89, 0x0b8d, 0x0b93, 0x0b99, 0x0b9b, 0x0bb7, 0x0bb9, 0x0bc3,
    0x0bcb, 0x0bcf, 0x0bdd, 0x0be1, 0x0be9, 0x0bf5, 0x0bfb, 0x0c07,
    0x0c0b, 0x0c11, 0x0c25, 0x0c2f, 0x0c31, 0x0c41, 0x0c5b, 0x0c5f,
    0x0c61, 0x0c6d, 0x0c73, 0x0c77, 0x0c83, 0x0c89, 0x0c91, 0x0c95,
    0x0c9d, 0x0cb3, 0x0cb5, 0x0cb9, 0x0cbb, 0x0cc7, 0x0ce3, 0x0ce5,
    0x0ceb, 0x0cf1, 0x0cf7, 0x0cfb, 0x0d01, 0x0d03, 0x0d0f, 0x0d13,
    0x0d1f, 0x0d21, 0x0d2b, 0x0d2d, 0x0d3d, 0x0d3f, 0x0d4f, 0x0d55,
    0x0d69, 0x0d79, 0x0d81, 0x0d85, 0x0d87, 0x0d8b, 0x0d8d, 0x0da3,
    0x0dab, 0x0db7, 0x0dbd, 0x0dc7, 0x0dc9, 0x0dcd, 0x0dd3, 0x0dd5,
    0x0ddb, 0x0de5, 0x0de7, 0x0df3, 0x0dfd, 0x0dff, 0x0e09, 0x0e17,
    0x0e1d, 0x0e21, 0x0e27, 0x0e2f, 0x0e35, 0x0e3b, 0x0e4b, 0x0e57,
    0x0e59, 0x0e5d, 0x0e6b, 0x0e71, 0x0e75, 0x0e7d, 0x0e87, 0x0e8f,
    0x0e95, 0x0e9b, 0x0eb1, 0x0eb7, 0x0eb9, 0x0ec3, 0x0ed1, 0x0ed5,
    0x0edb, 0x0eed, 0x0eef, 0x0ef9, 0x0f07, 0x0f0b, 0x0f0d, 0x0f17,
    0x0f25, 0x0f29, 0x0f31, 0x0f43, 0x0f47, 0x0f4d, 0x0f4f, 0x0f53,
    0x0f59, 0x0f5b, 0x0f67, 0x0f6b, 0x0f7f, 0x0f95, 0x0fa1, 0x0fa3,
    0x0fa7, 0x0fad, 0x0fb3, 0x0fb5, 0x0fbb, 0x0fd1, 0x0fd3, 0x0fd9,
    0x0fe9, 0x0fef, 0x0ffb, 0x0ffd, 0x0000
};

/* 564 primes < 2^12 factor 86.53% of all odd integers. */

/******************************************************************************/

static inline unsigned int sp_factor (uint32_t p[], uint32_t n)
{
    uint32_t sp = sp_lut[0], q;
    unsigned int np = 1;

    /* assert(n > 1 && n < (UINT32_C(1) << (24))); */

    for (unsigned int i = 1; (q = n / sp) >= sp; )
    {
        if (q * sp == n)
            np++, n = q, *p++ = sp;

        else if ((sp = sp_lut[i++]) == 0) /* EOT entry: */
            break;
    }

    *p++ = n;

    return np; /* the number of prime factors (with multiplicity). */
}

/******************************************************************************/

#endif /* SP_K12_H_ */
