#ifndef _BBOX_H_
#define _BBOX_H_

#include <algorithm>

#include "..\vec3.h"
#include "..\constant.h"

#include "ray.h"

namespace hstd {

namespace rt {
    
class BBox {
private:
public:
    Float3 pmin, pmax;

    BBox() {
        pmin = Float3(kINF, kINF, kINF);
        pmax = -1.0f * pmin;
    }

    BBox(const Float3 &p) : pmin(p), pmax(p) {}

    BBox(const Float3 &p1, const Float3 &p2) {
        pmin = Float3(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
        pmax = Float3(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
    }

    bool inside(const Float3 &pt) const {
        return (pmin.x <= pt.x && pt.x <= pmax.x &&
                pmin.y <= pt.y && pt.y <= pmax.y &&
                pmin.z <= pt.z && pt.z <= pmax.z);
    }

    Float3 &operator[](int i) {
        if (i == 0)
            return pmin;
        return pmax;
    }	
    const Float3 &operator[](int i) const {
        if (i == 0)
            return pmin;
        return pmax;
    }

    void expand(float delta) {
        const Float3 v(delta, delta, delta);
        pmin = pmin - v;
        pmax = pmax + v;
    }

    float surfaceArea() {
        const Float3 d = pmax - pmin;
        return 2.0f * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    float volume() {
        const Float3 d = pmax - pmin;
        return d.x * d.y * d.z;
    }

    enum LongestAxis {
        AxisX = 0,
        AxisY,
        AxisZ,
    };

    LongestAxis maximumExtent() const {
        const Float3 diag = pmax - pmin;
        if (diag.x > diag.y && diag.x > diag.z)
            return AxisX;
        else if (diag.y > diag.z)
            return AxisY;
        else
            return AxisZ;
    }

    // return true, if intersected
    bool checkIntersect(const Ray &ray, float *hitt0, float *hitt1) const {
        float t0 = 0.0, t1 = kINF;
        for (int i = 0; i < 3; ++i) {
            // Update interval for _i_th bounding box slab
            float invRayDir = 1.f / ray.dir[i];
            float tNear = (pmin[i] - ray.org[i]) * invRayDir;
            float tFar  = (pmax[i] - ray.org[i]) * invRayDir;

            // Update parametric interval from slab intersection $t$s
            if (tNear > tFar) std::swap(tNear, tFar);
            t0 = tNear > t0 ? tNear : t0;
            t1 = tFar  < t1 ? tFar  : t1;
            if (t0 > t1) return false;
        }
        if (hitt0) *hitt0 = t0;
        if (hitt1) *hitt1 = t1;
        return true;
    }

    friend BBox unionBBox(const BBox &b, const Float3 &p);
    friend BBox unionBBox(const BBox &b1, const BBox &b2);
};


inline BBox unionBBox(const BBox &b, const Float3 &p) {
    BBox ret = b;
    ret.pmin.x = std::min(b.pmin.x, p.x);
    ret.pmin.y = std::min(b.pmin.y, p.y);
    ret.pmin.z = std::min(b.pmin.z, p.z);
        
    ret.pmax.x = std::max(b.pmax.x, p.x);
    ret.pmax.y = std::max(b.pmax.y, p.y);
    ret.pmax.z = std::max(b.pmax.z, p.z);

    return ret;
}

inline BBox unionBBox(const BBox &b1, const BBox &b2) {
    BBox ret = b1;
    ret.pmin.x = std::min(b1.pmin.x, b2.pmin.x);
    ret.pmin.y = std::min(b1.pmin.y, b2.pmin.y);
    ret.pmin.z = std::min(b1.pmin.z, b2.pmin.z);

    ret.pmax.x = std::max(b1.pmax.x, b2.pmax.x);
    ret.pmax.y = std::max(b1.pmax.y, b2.pmax.y);
    ret.pmax.z = std::max(b1.pmax.z, b2.pmax.z);

    return ret;
}

inline std::ostream &operator<<(std::ostream &out, const BBox &b) {
    out << "{" <<  b.pmin << " - " << b.pmax << "}";
    return out;
}


} // namespace rt

} // namespace hstd

#endif