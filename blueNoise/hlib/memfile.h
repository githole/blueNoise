#ifndef _MEM_FILE_H_
#define _MEM_FILE_H_

#include <memory>
#include <algorithm>
#include <vector>
#include <iostream>

namespace hstd {

// メモリに全部展開しちゃうよ
class FileManager {
public:
	typedef std::shared_ptr<FILE> FilePtr;

private:
	FilePtr fp_;
	std::vector<unsigned char> buffer_;
	int now_index_;

	bool skip_if_return() {
		if (now_index_ >= buffer_.size())
			return false;

		if (now_index_ + 1 < buffer_.size()) {
			if (buffer_[now_index_] == '\r' && buffer_[now_index_ + 1] == '\n') {
				now_index_ += 2;
				return true;
			}
		}

		if (buffer_[now_index_] == '\r' || buffer_[now_index_] == '\n') {
			now_index_ += 1;
			return true;
		}
		
		return false;
	}
public:
	FileManager() {}
	virtual ~FileManager() {}

	std::vector<unsigned char>& buffer() { return buffer_; }

	static bool save(const char *filename, unsigned char *data, size_t num_element) {
		FilePtr fp = FilePtr(fopen(filename, "wb"), [](FILE *f){ if (f != NULL) fclose(f); });
		 
		if (fp == NULL) {
			std::cerr << "Load Error: " << filename << std::endl;
			return false;
		}

		size_t left_element = num_element;

		for (;left_element > 0;) {
			const size_t wrote_element = fwrite(data + (num_element - left_element), sizeof(unsigned char), left_element, fp.get());
			if (wrote_element <= 0)
				return false;
			left_element -= wrote_element;
		}

		return true;
	}

	bool load(const char *filename) {
		 fp_ = FilePtr(fopen(filename, "rb"), [](FILE *f){ if (f != NULL) fclose(f); });
		 
		if (fp_ == NULL) {
			std::cerr << "Load Error: " << filename << std::endl;
			return false;
		}
		
		// いまどきメモリたくさんあるっしょ
		// ってことで一括読み込み
		const long now_pos = ftell(fp_.get());
		fseek(fp_.get(), 0, SEEK_END);
		const long end_pos = ftell(fp_.get());
		fseek(fp_.get(), now_pos, SEEK_SET);

		const long total_size = end_pos - now_pos;
		buffer_.resize(total_size);

		
		const size_t ret_size = fread(&buffer_[0], sizeof(unsigned char), total_size, fp_.get());
		if (ret_size < total_size) {
			std::cerr << "Error : fread" << std::endl;
			return false;
		}

		now_index_ = 0;
		return true;
	}

	// 一行読み込み
	// 改行はなし
	bool gets(std::string* tmp) {
		if (now_index_ >= buffer_.size())
			return false;


		for (;;) {
			if (skip_if_return()) 
				break;

			if (now_index_ < buffer_.size())
				*tmp += buffer_[now_index_];

			++now_index_;
		}

		return true;
	}
};

};

#endif // _MEM_FILE_H_