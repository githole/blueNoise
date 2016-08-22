#include "hlib/image.h"
#include "hlib/random.h"


void Swap(hstd::Image& image, int px, int py, int qx, int qy)
{
    hstd::Color tmp = image.at(px, py);
    image.at(px, py) = image.at(qx, qy);
    image.at(qx, qy) = tmp;
}

const int W = 128;
const int H = 128;

const float sigma_i = 2.1f;
const float sigma_s = 1.0f;

float Comp(hstd::Image& image, int ix, int iy)
{
    float sum[H] = {};


    // small kernel approx.
#pragma omp parallel for schedule(dynamic, 1) num_threads(7)
    for (int oy = -7; oy <= 7; ++oy) {
        for (int ox = -7; ox <= 7; ++ox) {
            int sx = ix + ox;
            if (sx < 0)
                sx += W;
            if (sx >= W)
                sx -= W;

            int sy = iy + oy;
            if (sy < 0)
                sy += H;
            if (sy >= H)
                sy -= H;

            float dx = abs(ix - sx);
            if (dx > W / 2)
                dx = W - dx;

            float dy = abs(iy - sy);
            if (dy > H / 2)
                dy = H - dy;
            const float a =
                (dx * dx + dy * dy) / (sigma_i * sigma_i);

#ifdef OneDim
            const float b =
                sqrt(abs(image.at(ix, iy).x - image.at(sx, sy).x)) / (sigma_s * sigma_s);
#else
            const float da = abs(image.at(ix, iy)[0] - image.at(sx, sy)[0]);
            const float db = abs(image.at(ix, iy)[1] - image.at(sx, sy)[1]);

            const float b =
                sqrt(da * da + db * db) / (sigma_s * sigma_s);
#endif

            sum[sy] += exp(-a - b);
        }
    }

    /*
    // full size kernel
#pragma omp parallel for schedule(dynamic, 1) num_threads(7)
    for (int sy = 0; sy < H; ++sy) {
        for (int sx = 0; sx < W; ++sx) {

            float dx = abs(ix - sx);
            if (dx > W / 2)
                dx = W - dx;

            float dy = abs(iy - sy);
            if (dy > H / 2)
                dy = H - dy;
            const float a =
                (dx * dx + dy * dy) / (sigma_i * sigma_i);
            
            const float b =
                sqrt(abs(image.at(ix, iy).x - image.at(sx, sy).x)) / (sigma_s * sigma_s);

            sum[sy] += exp(-a - b);
        }
    }
    */

    float total = 0;
    for (int sy = 0; sy < H; ++sy)
        total += sum[sy];
    return total;
}

int main()
{
    hstd::Image image(W, H);

    hstd::Random random(0);

    // init
    for (int iy = 0; iy < H; ++iy) {
        for (int ix = 0; ix < W; ++ix) {
#ifdef OneDim
            const float r = random.next01();
            image.at(ix, iy) = hstd::Color(r, r, r);
#else
            const float r = random.next01();
            const float g = random.next01();
            image.at(ix, iy) = hstd::Color(r, g, 0);
#endif
        }
    }

    hstd::HDROperator::save("noise.hdr", &image);

    // initial energy
    float energy[H][W] = {};

    for (int iy = 0; iy < H; ++iy) {
        for (int ix = 0; ix < W; ++ix) {
            energy[iy][ix] = Comp(image, ix, iy);
        }
    }


    // random walk
    const int kMaxIteration = 1000000;

    for (int ite = 0; ite < kMaxIteration; ++ite) {
        float current_energy = 0;

        for (int iy = 0; iy < H; ++iy) {
            for (int ix = 0; ix < W; ++ix) {
                current_energy += energy[iy][ix];
            }
        }
        // printf("[%f]", current_energy);

        if (ite % 100000 == 0)
            printf("%d, %f\n", ite, current_energy);


        const int px = random.next() % W;
        const int py = random.next() % H;
        const int qx = random.next() % W;
        const int qy = random.next() % H;

        float next_energy = current_energy;
        next_energy -= energy[py][px];
        next_energy -= energy[qy][qx];

        // test swap
        Swap(image, px, py, qx, qy);

        const float e0 = Comp(image, px, py);
        const float e1 = Comp(image, qx, qy);

        next_energy += (e0 + e1);

        if (next_energy < current_energy) {
            energy[py][px] = e0;
            energy[qy][qx] = e1;
            continue;
        }

        // recover
        Swap(image, px, py, qx, qy);
    }

    hstd::HDROperator::save("result.hdr", &image);

    return 0;
}