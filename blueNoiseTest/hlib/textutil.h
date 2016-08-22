#ifndef _TEXTUTIL_H_
#define _TEXTUTIL_H_

#include <vector>
#include <string>

namespace hstd {

std::vector<std::string> split(const std::string& input, const char delimiter) {
	size_t now_index = 0;
	std::vector<std::string> ret;
	std::string tmp;

	for (;now_index < input.length(); ++now_index) {
		if (input[now_index] == delimiter) {
			if (tmp != "") {
				ret.push_back(tmp);
				tmp = "";
			}
			continue;
		}

		tmp += input[now_index];
	}
	if (tmp != "") {
		ret.push_back(tmp);
		tmp = "";
	}
	return ret;
}

};

#endif // _TEXTUTIL_H_
