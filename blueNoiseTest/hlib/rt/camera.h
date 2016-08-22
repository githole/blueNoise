#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "../vec3.h"
#include "ray.h"

namespace hstd {

namespace rt {

// ÉJÉÅÉâ
// TODO: ÉÅÉìÉoïœêîÇÃñºëOÇ™Ç§ÇÒÇ±
class Camera {
public:
	Float3 camera_direction_;
	Float3 camera_lookat_;
	Float3 camera_position_;

	Float3 screen_center_;
	Float3 screen_side_;
	Float3 screen_up_;

	float screen_width_;
	float screen_height_;

	float camera_near_, camera_far_;

	Camera(Float3& camera_position, Float3& camera_lookat, Float3& camera_up, float screen_width, float screen_height, float camera_near, float camera_far) : 
	camera_position_(camera_position),
	camera_lookat_(camera_lookat),
	screen_width_(screen_width),
	screen_height_(screen_height),
	camera_near_(camera_near),
	camera_far_(camera_far) {
		camera_direction_ = normalize(camera_lookat_ - camera_position_);
		
		screen_side_ = normalize(cross(camera_direction_, camera_up));
		screen_up_   = normalize(cross(camera_direction_, screen_side_));
		screen_center_ = camera_position_;
	}

	float screen_width() { return screen_width_; }
	float screen_height() { return screen_height_; }

	virtual rt::Ray getRay(const float u, const float v) const = 0;
	virtual Float3 world2screenUV(const Float3 &pos) = 0;
};

} // namespace rt

} // namespace hstd


#endif // _CAMERA_H_