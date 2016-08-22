#ifndef _IMAGE_H_
#define _IMAGE_H_

#include <memory>
#include <algorithm>
#include <vector>
#include "vec3.h"

#include "./ext/tiffio.h"

namespace hstd {

typedef Float3 Color;

class Image {
private:
	std::vector<Color> data_;
	unsigned int width_, height_;
public:
	explicit Image() :
	width_(0), height_(0), data_(0) {}

	explicit Image(unsigned int width, unsigned int height) :
	width_(width), height_(height), data_(width * height) {}

	void resize(unsigned int width, unsigned int height) {
		width_ = width;
		height_ = height;
		data_.resize(width_ * height_, Color(0, 0, 0));
	}

	unsigned int width() const { 
		return width_;
	}
	unsigned int height() const {
		return height_;
	}

	Color& sampleLoop(unsigned int x, unsigned int y) {
		unsigned int ix = x % width_;
		unsigned int iy = y % height_;

		return at(ix, iy);
	}

	inline unsigned int index(unsigned int x, unsigned int y) const {
		const unsigned int i = y * width_ + x;
#ifdef _DEBUG
		assert(i < data_.size());
#endif // _DEBUG
		return i;
	}

	Color& at(unsigned int x, unsigned int y) {
		return data_[index(x, y)];
	};
	
	const Color& at(unsigned int x, unsigned int y) const {
		return data_[index(x, y)];
	};
};


class TiffOperator {
private:

	template<typename T>
	static void read_(Image *image, uint32_t w, uint32_t h, void *raster) {
		T *arr = (T*)raster;
		T tmax =  std::numeric_limits<T>::max();
		for (int32_t y = 0; y < h; ++y) {
			for (int32_t x = 0; x < w; ++x) {
				Color col;
				// [0, 1]にスケーリング
				col.x = (float)arr[(y * w + x) * 3 + 0] / tmax;
				col.y = (float)arr[(y * w + x) * 3 + 1] / tmax;
				col.z = (float)arr[(y * w + x) * 3 + 2] / tmax;
				image->at(x, y) = col;
			}
		}
	}


public:
	static bool load(const char *filename, Image *image) {
		TIFF *tif = TIFFOpen(filename, "r");
		if (tif == NULL) {
			std::cout << "x" << std::endl;
			return false;
		}

		uint32 w, h;
		uint16 bps;
		size_t npixels;
		void *raster;
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
		npixels = w * h;

		uint32_t bytes_per_ch = bps / 8;

		image->resize(w, h);

		raster = _TIFFmalloc(3 * npixels * sizeof (bytes_per_ch));

		if (raster != NULL) {
		
			for (int32_t l = 0; l < h; ++l) {
				TIFFReadScanline(tif, (uint8_t*)raster + l * w * 3 * bytes_per_ch, l, 0);
			}

			if (bytes_per_ch == 1) {
				read_<uint8_t>(image, w, h, raster);
			} else if (bytes_per_ch == 2) {
				read_<uint16_t>(image, w, h, raster);
			}
		
			_TIFFfree(raster);
		}
		TIFFClose(tif);

		return true;
	}
	template<typename TypePerSample>
	static void save(const char *filename, const Image *image, const float min_value, const float max_value) {
		std::shared_ptr<TIFF> fp(TIFFOpen(filename, "w"), [](TIFF *f){ if (f != NULL) TIFFClose(f); });
		if (fp == NULL) {
			std::cerr << "Error: " << filename << std::endl;
			return;
		}

		TIFFSetField(fp.get(), TIFFTAG_IMAGEWIDTH, image->width());
		TIFFSetField(fp.get(), TIFFTAG_IMAGELENGTH, image->height());
		TIFFSetField(fp.get(), TIFFTAG_BITSPERSAMPLE, sizeof(TypePerSample) * 8);
		TIFFSetField(fp.get(), TIFFTAG_SAMPLESPERPIXEL, 3); // RGB
		TIFFSetField(fp.get(), TIFFTAG_ROWSPERSTRIP, image->height());
		TIFFSetField(fp.get(), TIFFTAG_COMPRESSION, COMPRESSION_NONE);
		TIFFSetField(fp.get(), TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
		TIFFSetField(fp.get(), TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
		TIFFSetField(fp.get(), TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(fp.get(), TIFFTAG_XRESOLUTION, 72.0);
		TIFFSetField(fp.get(), TIFFTAG_YRESOLUTION, 72.0);
		TIFFSetField(fp.get(), TIFFTAG_RESOLUTIONUNIT, RESUNIT_CENTIMETER);
		
		const TypePerSample target_min_value = std::numeric_limits<TypePerSample>::min();
		const TypePerSample target_max_value = std::numeric_limits<TypePerSample>::max();

		std::vector<TypePerSample> data(image->width() * image->height() * 3);

		auto mp = [=](float x) {
			const float y = (x - min_value)  / (max_value - min_value);
			return (TypePerSample)(y * (target_max_value - target_min_value) + target_min_value);
		};

		for (int iy = 0; iy < image->height(); ++iy) {
			for (int ix = 0; ix < image->width(); ++ix) {
				const TypePerSample r = mp(image->at(ix, iy).x);
				const TypePerSample g = mp(image->at(ix, iy).y);
				const TypePerSample b = mp(image->at(ix, iy).z);
				data[(iy * image->width() + ix) * 3 + 0] = r;
				data[(iy * image->width() + ix) * 3 + 1] = g;
				data[(iy * image->width() + ix) * 3 + 2] = b;
			}
		}
		
		TIFFWriteEncodedStrip(fp.get(), 0, &data[0], sizeof(TypePerSample) * image->width() * image->height() * 3);
	}
};


class HDROperator {
private:
	struct HDRPixel {
		unsigned char r, g, b, e;
		HDRPixel(const unsigned char r = 0, const unsigned char g = 0, const unsigned char b = 0, const unsigned char e = 0) :
		r(r), g(g), b(b), e(e) {};
		unsigned char get(int idx) const {
			switch (idx) {
			case 0: return r;
			case 1: return g;
			case 2: return b;
			case 3: return e;
			} return 0;
		}

		HDRPixel(const Color& color) {
			double d = std::max(color.x, std::max(color.y, color.z));
			if (d <= 1e-32) {
				r = g = b = e = 0;
				return;
			}
			int ie;
			const double m = frexp(d, &ie); // d = m * 2^e
			d = m * 256.0 / d;

			r = (unsigned char)(color.x * d);
			g = (unsigned char)(color.y * d);
			b = (unsigned char)(color.z * d);
			e = (unsigned char)(ie + 128);
		}
	};
public:

	static bool load(const char *filename, Image *image) {
		std::shared_ptr<FILE> fp(fopen(filename, "rb"), [](FILE *f){ if (f != NULL) fclose(f); });

		const int BufSize = 4096;
		char buf[BufSize];
		if (fp == NULL) {
			std::cerr << "HDR Error: " << filename << std::endl;
			return false;
		}

		bool valid = false;
		enum FileType {
			NONE,
			RLE_RGBE_32,
		} type;

		type = NONE;
		float exposure = 1.0;

		// ヘッダ読み込み
		for (;;) {
			fgets(buf, BufSize, fp.get());
			// std::cout << buf;
			if (buf[0] == '#') {
				if (strcmp(buf, "#?RADIANCE\n") == 0)
					valid = true;
			} else {
				if (strstr(buf, "FORMAT=") == buf) {
					char buf2[BufSize];
					sscanf(buf, "FORMAT=%s", buf2);
					if (strcmp(buf2, "32-bit_rle_rgbe") == 0)
						type = RLE_RGBE_32;

				} else if (strstr(buf, "EXPOSURE=") == buf) {
					sscanf(buf, "FORMAT=%f", &exposure);
				}
			} 
			
			if (buf[0] == '\n')
				break;
		}

		if (!valid) {
			std::cerr << "Invalid HDR File: " << filename << std::endl;
			return false;
		}

		// サイズ読み込み
		char buf2[BufSize], buf3[BufSize];
		int width, height;
		fgets(buf, BufSize, fp.get());
		sscanf(buf, "%s %d %s %d", buf2, &height, buf3, &width);

		// -Y @ +X @ しかあつかわないことにする　
		if (strcmp(buf2, "-Y") != 0 && strcmp(buf3, "+X") != 0) {
			std::cerr << "Invalid HDR File: " << filename << " " << buf2 << " " << buf3 << std::endl;
			return false;
		}

		std::vector<unsigned char> tmp_data(width * height * 4);


		long now_pos = ftell(fp.get());
		fseek(fp.get(), 0, SEEK_END);
		long end_pos = ftell(fp.get());
		fseek(fp.get(), now_pos, SEEK_SET);

		const long rest_size = end_pos - now_pos;
		std::vector<unsigned char> buffer(rest_size);

		const size_t ret_size = fread(&buffer[0], sizeof(unsigned char), rest_size, fp.get());
		if (ret_size < rest_size) {
			std::cerr << "Error: fread" << std::endl;
			return false;
		}
		
		int index = 0;
		int nowy = 0;
		for (;index < rest_size;) {
			const int now = buffer[index++];
			if (now == EOF)
				break;
			const int now2 = buffer[index++];

			if (now != 0x02 || now2 != 0x02) {
				break;
			}

			const int A = buffer[index++];
			const int B = buffer[index++];
			const int width = (A << 8) + B;

			int nowx = 0;
			int nowvalue = 0;
			for (;;) {
				if (nowx >= width) {
					nowvalue ++;
					nowx = 0;
					if (nowvalue == 4)
						break;
				}

				const int info = buffer[index++];
				if (info <= 128) {
					for (int i = 0; i < info; ++i) {
						const int data = buffer[index++];
						tmp_data[(nowy * width + nowx) * 4 + nowvalue] = data;
						nowx ++;
					}
				} else {
					const int num = info - 128;
					const int data = buffer[index++];
					for (int i = 0; i < num; ++i) {
						tmp_data[(nowy * width + nowx) * 4 + nowvalue] = data;
						nowx ++;
					}
				}
			}

			nowy ++;
		}

		image->resize(width, height);

		// 展開
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				const int e = tmp_data[(y * width + x) * 4 + 3];
				image->at(x, y).x = tmp_data[(y * width + x) * 4 + 0] * pow(2, e - 128.0f) / 256.0f;
				image->at(x, y).y = tmp_data[(y * width + x) * 4 + 1] * pow(2, e - 128.0f) / 256.0f;
				image->at(x, y).z = tmp_data[(y * width + x) * 4 + 2] * pow(2, e - 128.0f) / 256.0f;
			}
		}

		return true;
	}

	static void save(const char *filename, Image *image, bool absolute = false) {
		std::shared_ptr<FILE> fp(fopen(filename, "wb"), [](FILE *f){ if (f != NULL) fclose(f); });

		if (fp == NULL) {
			std::cerr << "Error: " << filename << std::endl;
			return;
		}
		// .hdrフォーマットに従ってデータを書きだす
		// ヘッダ
		unsigned char ret = 0x0a;
		fprintf(fp.get(), "#?RADIANCE%c", (unsigned char)ret);
		fprintf(fp.get(), "# Made with 100%% pure HDR Shop%c", ret);
		fprintf(fp.get(), "FORMAT=32-bit_rle_rgbe%c", ret);
		fprintf(fp.get(), "EXPOSURE=1.0000000000000%c%c", ret, ret);

		// 輝度値書き出し
		const unsigned int width = image->width();
		const unsigned int height = image->height();

		fprintf(fp.get(), "-Y %d +X %d%c", height, width, ret);

		std::vector<unsigned char> buffer;
		
		for (int i = 0; i < height; ++ i) {
			std::vector<HDRPixel> line;
			for (int j = 0; j < width; ++ j) {
				Color color = image->at(j, i);
				if (absolute) {
					color.x = abs(color.x);
					color.y = abs(color.y);
					color.z = abs(color.z);
				}
				HDRPixel p(color);
				line.push_back(p);
			}
			buffer.push_back(0x02);
			buffer.push_back(0x02);
			buffer.push_back((width >> 8) & 0xFF);
			buffer.push_back(width & 0xFF);
			for (int i = 0; i < 4; i ++) {
				for (int cursor = 0; cursor < width;) {
					const int cursor_move = std::min((unsigned int)127, width - cursor);
					buffer.push_back(cursor_move);
					for (int j = cursor; j < cursor + cursor_move; j ++) {
						buffer.push_back(line[j].get(i));
					}
					cursor += cursor_move;
				}
			}
		}
		fwrite(&buffer[0], sizeof(unsigned char), buffer.size(), fp.get());
	}
};




};

#endif // _IMAGE_H_