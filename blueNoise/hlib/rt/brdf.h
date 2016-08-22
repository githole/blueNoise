#ifndef _BRDF_H_
#define _BRDF_H_

#include "../image.h"
#include "../random.h"
#include "../constant.h"
#include "../sampling.h"
#include "../hmath.h"

namespace hstd {

namespace rt {
	

class BRDF {
private:
public:
	virtual Color reflectance() const = 0;
	virtual Color eval(const Float3& in, const Float3& normal, const Float3& out) const = 0;
	virtual Float3 sample(Random& random, const Float3& in, const Float3& normal, float* pdf) const = 0;
	virtual float eval_pdf(const Float3& in, const Float3& normal, const Float3& out) const = 0;
};

class LambertianBRDF : public BRDF {
private:
	Color reflectance_;
public:
	virtual Color reflectance() const {
		return reflectance_;
	}

	LambertianBRDF(const Color& reflectance) : reflectance_(reflectance) {}
	virtual Color eval(const Float3& in, const Float3& normal, const Float3& out) const {
		return reflectance_ / kPI;
	}

	virtual float eval_pdf(const Float3& in, const Float3& normal, const Float3& out) const {
		return dot(normal, out) / kPI;
	}

	virtual Float3 sample(Random& random, const Float3& in, const Float3& normal, float* pdf) const {
		Float3 binormal, tangent, now_normal = normal;
			
		/*
		if (dot(in, normal) > 0)
			now_normal = -now_normal;
		*/
							
		createOrthoNormalBasis(now_normal, &tangent, &binormal);
		const Float3 dir = Sampling::cosineWeightedHemisphereSurface(random, now_normal, tangent, binormal);
		// cosΘ/pi
		if (pdf != NULL) {
			*pdf = eval_pdf(in, now_normal, dir);
		}
		return dir;
	}
};


class PhongBRDF : public BRDF {
private:
	Color reflectance_;
	float n_;
public:
	PhongBRDF(const Color& reflectance, const float n) : reflectance_(reflectance), n_(n) {}
	
	virtual Color reflectance() const {
		return reflectance_;
	}

	virtual Color eval(const Float3& in, const Float3& normal, const Float3& out) const {
		Float3 reflection_dir = reflect(in, normal);
		float cosa = dot(reflection_dir, out);
		if (cosa < 0)
			cosa = 0.0f;
		return reflectance_ * (n_ + 2.0f) / (2.0f * kPI) * pow(cosa, n_);
	}

	virtual float eval_pdf(const Float3& in, const Float3& normal, const Float3& out) const {
		Float3 reflection_dir = reflect(in, normal);
		float cosa = dot(reflection_dir, out);
		if (cosa < 0)
			cosa = 0.0f;
		return (n_ + 1.0f) / (2.0f * kPI) * pow(cosa, n_);
	}

	virtual Float3 sample(Random& random, const Float3& in, const Float3& normal, float* pdf) const {
		Float3 dir;
		Float3 reflection_dir = reflect(in, normal);
		Float3 binormal, tangent;
		createOrthoNormalBasis(reflection_dir, &tangent, &binormal);

		float u1 = random.next01();
		float u2 = random.next01();

		float theta = acos(pow(u1, 1 / (n_ + 1)));
		float phi = u2 * 2.0f * kPI;
							
		dir = tangent * sin(theta) * cos(phi) + reflection_dir * cos(theta) + binormal *sin(theta) * sin(phi);
		
		if (pdf != NULL) {
			*pdf = eval_pdf(in, normal, dir);
		}

		return dir;
	}
};



} // namespace rt

} // namespace hstd

#endif // _BRDF_H_