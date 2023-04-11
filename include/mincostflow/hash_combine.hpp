#pragma once

// Copyright 2005-2014 Daniel James.
// Copyright 2021 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Based on Peter Dimov's proposal
// http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2005/n1756.pdf
// issue 6.18.
//
// This also contains public domain code from MurmurHash. From the
// MurmurHash header:

// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

#include <cstdint>

namespace ln
{
    uint64_t hash_combine(uint64_t h, uint64_t k)
    {
        const uint64_t m = (uint64_t(0xc6a4a793) << 32) + 0x5bd1e995;
        const int r = 47;
        
        k*=m;
        k^= k>>r;
        k*=m;
        
        h^=k;
        h*=m;
        
        h += 0xe6546b64;
        return h;
    }
}
