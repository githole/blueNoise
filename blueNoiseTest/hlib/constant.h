#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#include <limits>

namespace hstd {


template<typename T>
inline T PI() {
	return 3.1415926535897931;
};

template<>
inline float PI() {
	return 3.1415927f;
};

template<>
inline double PI() {
	return 3.1415926535897931;
};

static const float kPI = PI<float>();

static const float kINF = std::numeric_limits<float>::infinity();

} // namespace hstd

#endif // _CONSTANT_H_
