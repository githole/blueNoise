#ifndef _IBL_H_
#define _IBL_H_

#include "../image.h"
#include "../random.h"
#include "../constant.h"
#include "../sampling.h"
#include "../hmath.h"
#include "brdf.h"

namespace hstd {
namespace rt {

class ImageBasedLight {
    
private:
    Image image_;

    Image importance_map_;
    
    Image warped_importance_map_;

    std::vector<float> volume_cdf;
    std::vector<float> volume_pdf;

    float luminance(const Color& color) const {
        return dot(Float3(0.2126f, 0.7152f, 0.0722), color);
    }

    static const int kWarpedImportanceMapSize = 64;
public:
    struct Sample {
        Float3 dir;
        float pdf;
    };

    /*
    ImageBasedLight(const char* filename) {
        HDROperator::load(filename, &image_);
    }*/

    bool load(const char *filename) {
        return HDROperator::load(filename, &image_);;
    }

    Color sample_from_direction(const Float3& dir) const {
        float theta, phi;
        directionToPolarCoordinate(dir, &theta, &phi);
                        
        const float iblv = theta / kPI;
        const float iblu = phi / (2.0f * kPI);
                        
        const float iblx = clamp(float(iblu * image_.width()), (float)0, image_.width() - 1.0f);
        const float ibly = clamp(float(iblv * image_.height()), (float)0, image_.height() - 1.0f);
                        
        return image_.at(iblx, ibly);
    }

    std::vector<Sample> create_sample_from_importance_map(const int num_sample, Random& random, const Float3& in, const Float3& normal,  BRDF& brdf) {
        Image now_importance_map = warped_importance_map_;

        const int width = now_importance_map.width();
        const int height = now_importance_map.height();
        
        Float3 now_normal = normal;
        if (dot(in, now_normal) > 0)
            now_normal = -now_normal;

        float total = 0;
        std::vector<float> pdf(width * height);
        std::vector<float> cdf(width * height);

        const float ox = random.next01();
        const float oy = random.next01();

        for (int iy = 0; iy < height; ++iy) {
            for (int ix = 0; ix < width; ++ix) {

                float value = 0;
                const int SuperSample = 1;
                for (int sy = 0; sy < SuperSample; ++sy) {
                    const float u = (float)(ix + ox) / width;
                    const float v = (float)(iy + oy) / height;

                    const float phi = u * 2.0f * kPI;
                    const float y = (1.0f - v) * 2.0f - 1.0f;

                    const float sqrt_y = sqrt(1 - y * y);

                    const Float3 dir = Float3(sqrt_y * cos(phi), y, sqrt_y * sin(phi));

                    const float cdot = clamp(dot(now_normal, dir), 0.0f, 1.0f);
        
                    value += cdot * luminance(brdf.eval(in, now_normal, dir)) / SuperSample;
                }

                now_importance_map.at(ix, iy) *= value;
//				now_importance_map.at(ix, iy) = Color(1, 1, 1);

                total += now_importance_map.at(ix, iy).x;
            }
        }

        // 手法A
        // CDFを作る
        pdf[0] = now_importance_map.at(0, 0).x / total;
        cdf[0] = pdf[0];
        for (int i = 1; i < width * height; ++i) {
            const int x = i % width;
            const int y = i / width;
            pdf[i] = now_importance_map.at(x, y).x / total;
            cdf[i] = pdf[i] + cdf[i - 1];
        }
        {
            std::vector<Sample> samples(num_sample);
            for (int i = 0; i < num_sample; ++i) {
                const float randomnumber = random.next01();
                std::vector<float>::iterator ite = lower_bound(cdf.begin(), cdf.end(), randomnumber);
                const int index = ite - cdf.begin();

                // 方向決定
                const int ix = index % width;
                const int iy = index / width;
                const float u = (float)(ix + random.next01()) / width;
                const float v = (float)(iy + random.next01()) / height;
                
                const float phi = u * 2.0f * kPI;
                const float y = (1.0f - v) * 2.0f - 1.0f;
                const Float3 out = Float3(sqrt(1 - y*y) * cos(phi), y, sqrt(1 - y*y) * sin(phi));

                //samples[i].dir = Sampling::uniformHemisphereSurface(random, normal, tangent, binormal);
                //samples[i].pdf = 1.0f / (2.0f * kPI);
                
                samples[i].dir = out;
                samples[i].pdf = pdf[index] * width * height * (1.0f / (4.0f * kPI));
                
            }

            return samples;
        }
    }

    
    Sample create_sample_volume(Random& random, const Float3& in) {
        Image now_importance_map = warped_importance_map_;

        const int width = now_importance_map.width();
        const int height = now_importance_map.height();
        
        Sample sample;

        const float randomnumber = random.next01();
        std::vector<float>::iterator ite = lower_bound(volume_cdf.begin(), volume_cdf.end(), randomnumber);
        const int index = clamp((int)(ite - volume_cdf.begin()), (int)0, (int)(volume_cdf.size() - 1));

        // 方向決定
        const int ix = index % width;
        const int iy = index / width;
        const float u = (float)(ix + random.next01()) / width;
        const float v = (float)(iy + random.next01()) / height;
                
        const float phi = u * 2.0f * kPI;
        const float y = (1.0f - v) * 2.0f - 1.0f;
        const Float3 out = Float3(sqrt(1 - y*y) * cos(phi), y, sqrt(1 - y*y) * sin(phi));

        /*
        if (index >= volume_pdf.size())
            std::cout << "*";
        */

        sample.dir = out;
        sample.pdf = volume_pdf[index] * width * height * (1.0f / (4.0f * kPI));

        return sample;
    }


    
    std::vector<Sample> create_sample_only_from_importance_map(const int num_sample, Random& random, const Float3& in, const Float3& normal, std::vector<float>& pdf) {
        Image now_importance_map = warped_importance_map_;

        const int width = now_importance_map.width();
        const int height = now_importance_map.height();
        
        Float3 now_normal = normal;
        if (dot(in, now_normal) > 0)
            now_normal = -now_normal;

        float total = 0;
        pdf.resize(width * height);
        std::vector<float> cdf(width * height);

        for (int iy = 0; iy < height; ++iy) {
            for (int ix = 0; ix < width; ++ix) {

                float value = 0;
                const int SuperSample = 1;
                for (int sy = 0; sy < SuperSample; ++sy) {
                    const float u = (float)(ix + random.next01()) / width;
                    const float v = (float)(iy + random.next01()) / height;

                    const float phi = u * 2.0f * kPI;
                    const float y = (1.0f - v) * 2.0f - 1.0f;

                    const Float3 dir = Float3(sqrt(1 - y*y) * cos(phi), y, sqrt(1 - y*y) * sin(phi));

                    const float cdot = clamp(dot(now_normal, dir), 0.0f, 1.0f);

                    value += cdot / SuperSample;
                }
                now_importance_map.at(ix, iy) *= value;
                total += now_importance_map.at(ix, iy).x;
            }
        }

        // 手法A
        // CDFを作る
        pdf[0] = now_importance_map.at(0, 0).x / total;
        cdf[0] = pdf[0];
        for (int i = 1; i < width * height; ++i) {
            const int x = i % width;
            const int y = i / width;
            pdf[i] = now_importance_map.at(x, y).x / total;
            cdf[i] = pdf[i] + cdf[i - 1];
        }
        {
            std::vector<Sample> samples(num_sample);
            for (int i = 0; i < num_sample; ++i) {
                const float randomnumber = random.next01();
                std::vector<float>::iterator ite = lower_bound(cdf.begin(), cdf.end(), randomnumber);
                const int index = ite - cdf.begin();

                // 方向決定
                const int ix = index % width;
                const int iy = index / width;
                const float u = (float)(ix + random.next01()) / width;
                const float v = (float)(iy + random.next01()) / height;
                
                const float phi = u * 2.0f * kPI;
                const float y = (1.0f - v) * 2.0f - 1.0f;
                const Float3 out = Float3(sqrt(1 - y*y) * cos(phi), y, sqrt(1 - y*y) * sin(phi));

                samples[i].dir = out;
                samples[i].pdf = pdf[index] * width * height * (1.0f / (4.0f * kPI));
            }
            

            // 確率密度関数マップを立体角測度に変換しとく（MIS用）
            for (int i = 1; i < width * height; ++i) {
                pdf[i] *= width * height * (1.0f / (4.0f * kPI));
            }

            return samples;
        }

    }

    // const float phi = u * 2.0f * kPI;
    // const float y = (1.0f - v) * 2.0f - 1.0f;
    // 
    // const Float3 dir = Float3(sqrt(1 - y*y) * cos(phi), y, sqrt(1 - y*y) * sin(phi));
    void index_from_direction(const Float3& dir, int *index) {
        float theta, phi;
        directionToPolarCoordinate(dir, &theta, &phi);

        const float v = 1.0f - (dir.y + 1.0f) / 2.0f;
        const float u = phi / (2.0f * kPI);

        const int x = clamp(int(u * importance_map_.width()), 0, int(importance_map_.width())-1);
        const int y = clamp(int(v * importance_map_.height()), 0, int(importance_map_.height())-1);

        *index = y * importance_map_.width() + x;
    }

    Float3 estimated_sundir;
    void create_importance_map(const int input_width, const int input_height) {
        importance_map_.resize(kWarpedImportanceMapSize, kWarpedImportanceMapSize);
        warped_importance_map_.resize(kWarpedImportanceMapSize, kWarpedImportanceMapSize);

        const int width = importance_map_.width();
        const int height = importance_map_.height();

        float total = 0;
        for (int iy = 0; iy < height; ++iy) {
            for (int ix = 0; ix < width; ++ix) {
                const float u = (float)ix / width;
                const float v = (float)iy / height;
                const float next_u = (float)(ix + 1) / width;
                const float next_v = (float)(iy + 1) / height;
                
                const int begin_x = u * image_.width();
                const int end_x = clamp(int(next_u * image_.width()), (int)0, (int)image_.width());
                const int begin_y = v * image_.height();
                const int end_y = clamp(int (next_v * image_.height()), (int)0, (int)image_.height());

                Color accum;
                const int area = (end_y - begin_y) * (end_x - begin_x);
                for (int ty = begin_y; ty < end_y; ++ty) {
                    for (int tx = begin_x; tx < end_x; ++tx) {
                        accum += image_.at(tx, ty);
                    }
                }

                importance_map_.at(ix, iy) = Float3(1, 1, 1) * luminance(accum) / area;
                total += importance_map_.at(ix, iy).x;
            }
        }
        // 正規化
        for (int iy = 0; iy < height; ++iy) {
            for (int ix = 0; ix < width; ++ix) {
                importance_map_.at(ix, iy) = importance_map_.at(ix, iy) / total;
            }
        }

        total = 0;
        // Warped
        // [0, 1]^2みたいな平面と球の間のマッピングをする
        // [-pi, pi]x[0, 2pi]みたいな平面から一様にサンプリングしても、球面上で一様なサンプリングにはならない。
        // [0, 1]^2から一様にサンプリングしたとき、球面上で一様にサンプリングするためにはワーピングしないといけない。
        const int SuperX = image_.width() / width;
        const int SuperY = image_.height() / height;
        for (int iy = 0; iy < height; ++iy) {
            for (int ix = 0; ix < width; ++ix) {
                Color accum;
                for (int sy = 0; sy < SuperY; ++sy) {
                    for (int sx = 0; sx < SuperX; ++sx) {
                        const float u = (float)(ix + (float)sx / SuperX) / width;
                        const float v = (float)(iy + (float)sy / SuperY) / height;

                        const float phi = u * 2.0f * kPI;
                        const float y = (1.0f - v) * 2.0f - 1.0f;

                        const Float3 dir = Float3(sqrt(1 - y*y) * cos(phi), y, sqrt(1 - y*y) * sin(phi));
                        accum += sample_from_direction(dir) / (SuperX * SuperY);
                    }
                }

                warped_importance_map_.at(ix, iy) = Float3(1, 1, 1) * luminance(accum);
                total += warped_importance_map_.at(ix, iy).x;
            }
        }

        // 正規化
        for (int iy = 0; iy < height; ++iy) {
            for (int ix = 0; ix < width; ++ix) {
                warped_importance_map_.at(ix, iy) = warped_importance_map_.at(ix, iy) / total;
            }
        }
        /*
        HDROperator::save("warped_importance.hdr", &warped_importance_map_);
        HDROperator::save("importance.hdr", &importance_map_);
        */

        // ボリュームレンダリング用
        
        volume_cdf.resize(width * height);
        volume_pdf.resize(width * height);

        // 手法A
        // CDFを作る
        volume_pdf[0] = warped_importance_map_.at(0, 0).x;
        volume_cdf[0] = volume_pdf[0];
        for (int i = 1; i < width * height; ++i) {
            const int x = i % width;
            const int y = i / width;
            volume_pdf[i] = warped_importance_map_.at(x, y).x;
            volume_cdf[i] = volume_pdf[i] + volume_cdf[i - 1];
        }

        {
            // 太陽の方向推定
            float vMax = -1;
            int index = -1;
            for (int i = 0; i < width * height; ++i) {
                if (vMax < volume_pdf[i]) {
                    vMax = volume_pdf[i];
                    index = i;
                }
            }
            // 方向決定
            const int ix = index % width;
            const int iy = index / width;
            const float u = (float)(ix + 0.5f) / width;
            const float v = (float)(iy + 0.5f) / height;
                
            const float phi = u * 2.0f * kPI;
            const float y = (1.0f - v) * 2.0f - 1.0f;
            estimated_sundir = Float3(sqrt(1 - y*y) * cos(phi), y, sqrt(1 - y*y) * sin(phi));
        }
    }
};

} // namespace rt
} // namespace hstd

#endif // _IBL_H_