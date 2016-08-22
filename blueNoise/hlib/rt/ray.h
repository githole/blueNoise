#ifndef _RAY_H_
#define _RAY_H_

#include "../vec3.h"

namespace hstd {

namespace rt {

struct Ray {
	Float3 org, dir;
	Ray(const Float3 &org, const Float3 &dir) : org(org), dir(dir) {}
	Ray() {}
};

}

}


#endif // _RAY_H_
