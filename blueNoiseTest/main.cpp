#include "hlib/image.h"
#include "hlib/rt/ray.h"
#include "hlib/rt/sphere.h"
#include "hlib/random.h"
#include "hlib/sampling.h"

hstd::Float3 next_dir(float u0, float u1, const hstd::Float3& normal, const hstd::Float3& tangent, const hstd::Float3& binormal)
{
    const float phi = u0 * 2.0f * hstd::kPI;
    const float r2 = u1, r2s = sqrt(r2);
          
    const float tx = r2s * cos(phi);
    const float ty = r2s * sin(phi);
    const float tz = sqrt(1.0f - r2);

    return tz * normal + tx * tangent + ty * binormal;
}

void offset(float& u0, float& u1, int x, int y, const hstd::Image& offsetImage)
{
    const int sx = x % offsetImage.width();
    const int sy = y % offsetImage.height();
    const auto& s = offsetImage.at(sx, sy);

    u0 += s[0];
    u1 += s[1];

    if (u0 >= 1)
        u0 -= 1;
    if (u1 >= 1)
        u1 -= 1;
}

//#define AO
#define MB
//#define BlueNoise

int main()
{
    const int W = 512;
    const int H = 512;
    hstd::Image image(W, H);

    // camera
    const hstd::Float3 org(0, 0.8f, -3);
    const hstd::Float3 dir = hstd::normalize(hstd::Float3(0, -0.2f, 1));
    const hstd::Float3 up(0, 1, 0);
    const float dist = 0.2f;
    const hstd::Float3 center = org + dist * dir;
    const hstd::Float3 screen_x = hstd::normalize(hstd::cross(dir, up));
    const hstd::Float3 screen_y = hstd::normalize(hstd::cross(screen_x, dir));
    const float screen_w = 0.1f;
    const float screen_h = 0.1f;

    // scene
    const int kSceneSize = 3;
    hstd::rt::Sphere spheres[kSceneSize] = {
        hstd::rt::Sphere(0.3f, hstd::Float3(-0.5, 0.3f, 0)),
        hstd::rt::Sphere(0.5f, hstd::Float3(0.3, 0.5f, 0)),
        hstd::rt::Sphere (10000.0f, hstd::Float3(0, -10000, 0)),
    };

    hstd::Random random(0);

    const int kSampleDiv = 1;
    const int kSample = kSampleDiv * kSampleDiv;
    const float kMin = 1e-3f;

    // generate samples
#ifdef AO
    float u0_seq[kSample];
    float u1_seq[kSample];

    for (int iy = 0; iy < kSampleDiv; ++iy) {
        for (int ix = 0; ix < kSampleDiv; ++ix) {
            const int index = ix + iy * kSampleDiv;
            u0_seq[index] = (ix + 0.5f) / kSampleDiv;
            u1_seq[index] = (iy + 0.5f) / kSampleDiv;
        }
    }
#else
    float u0_seq[kSample];
    for (int ix = 0; ix < kSample; ++ix) {
        u0_seq[ix] = (ix + 0.5f) / kSample;
    }
#endif

    // load offset
    hstd::Image offsetImage;
#ifdef BlueNoise
#ifdef AO
    hstd::HDROperator::load("bluenoise_2d.hdr", &offsetImage);
#else
    hstd::HDROperator::load("bluenoise_1d.hdr", &offsetImage);
#endif
#else
    hstd::HDROperator::load("noise.hdr", &offsetImage);
#endif

    const hstd::Float3 light(0, 3, 0);

    // check intersection
    auto check = [&spheres, &kSceneSize](const hstd::rt::Ray& ray, int& id, float& len) {
        for (int s = 0; s < kSceneSize; ++s) {
            hstd::rt::Hitpoint hp;
            auto& sphere = spheres[s];
            if (sphere.intersect(ray, &hp)) {
                if (len > hp.distance) {
                    len = hp.distance;
                    id = s;
                }
            }
        }
    };

    // rendering
    for (int iy = 0; iy < H; ++iy) {
        for (int ix = 0; ix < W; ++ix) {

            const float u = ((ix + 0.5f) / W) * 2 - 1;
            const float v = 1 - ((iy + 0.5f) / H) * 2;

            const hstd::Float3 screen_pos =
                center + u * screen_w * screen_x + v * screen_h * screen_y;
            const hstd::rt::Ray ray(org, hstd::normalize(screen_pos - org));


#ifdef MB
            for (int sample = 0; sample < kSample; ++sample) {
                float u0 = u0_seq[sample];
                float u1 = 0;
                offset(u0, u1, ix, iy, offsetImage);
                spheres[1].position.y = 0.5f + u0 * 0.3f;
                bool hit = false;
                float d = 1e+6f;
                int id = -1;

                check(ray, id, d);
                if (id >= 0) {
                    const hstd::Float3 pos = ray.org + d * ray.dir;
                    const hstd::Float3 normal = hstd::normalize(pos - spheres[id].position);

                    // shadow
                    const hstd::rt::Ray shadow_ray(pos + kMin * normal, hstd::normalize(light - pos));
                    float sd = 1e+6f;
                    int sid = -1;
                    check(shadow_ray, sid, sd);
                    if (sid == -1) {
                        // lighting
                        image.at(ix, iy) += std::max(0.0f, hstd::dot(normal, shadow_ray.dir)) * hstd::Float3(1, 1, 1);
                    }
                }
            }
            image.at(ix, iy) = image.at(ix, iy)/ kSample;
#endif

#ifdef AO

            bool hit = false;
            float d = 1e+6f;
            int id = -1;

            check(ray, id, d);

            if (id >= 0) {
                const hstd::Float3 pos = ray.org + d * ray.dir;
                const hstd::Float3 normal = hstd::normalize(pos - spheres[id].position);

                hstd::Float3 tangent, binormal;
                hstd::createOrthoNormalBasis(normal, &tangent, &binormal);

                // ao
                int numhit = 0;
                for (int sample = 0; sample < kSample; ++sample) {
                    float u0 = u0_seq[sample];
                    float u1 = u1_seq[sample];
                    offset(u0, u1, ix, iy, offsetImage);
                    const hstd::Float3 next = next_dir(u0, u1, normal, tangent, binormal);
                    const hstd::rt::Ray ao_ray(pos + kMin * normal, next);

                    float sd = 1e+6f;
                    int sid = -1;
                    check(ao_ray, sid, sd);

                    if (sid >= 0) {
                        ++numhit;
                    }
                }

                const float ao = 1 - (float)numhit / (kSample);
                image.at(ix, iy) = ao * hstd::Float3(1, 1, 1);
            }
#endif
        }
    }

    hstd::HDROperator::save("result.hdr", &image);

    return 0;
}