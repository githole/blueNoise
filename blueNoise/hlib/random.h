#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <limits>
#include <cstdint>

namespace hstd {

// Xor-Shiftによる乱数ジェネレータ
class Random {
private:
    uint64_t x;
public:
    uint64_t next(void) { 
        x ^= x >> 12; // a
        x ^= x << 25; // b
        x ^= x >> 27; // c
        return x * 2685821657736338717LL;
    }

    float next01(void) {
        return (float)((double)next() / std::numeric_limits<uint64_t>::max());
    }

    // [0, 1)
    float nextc01o(void) {
        return (float)((uint32_t)(next() >> 1) / 4294967296.0);
    }

    // [min_value, max_value]
    float next(float min_value, float max_value) {
        const double inv = (max_value - min_value);
        return (float)(((double)next() * (inv / std::numeric_limits<uint64_t>::max())) + min_value);
    }

    Random(uint64_t initial_seed) {
        if (initial_seed == 0)
            x = 0xDEADBEEFDEADBEEF; // xorshift64*のseedは非ゼロでないといけない。
        else
            x = initial_seed;
    }
};



};



#endif
