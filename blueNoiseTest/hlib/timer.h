#ifndef _TIMER_H_
#define _TIMER_H_


namespace hstd {
	
#ifdef _WIN32

#define NOMINMAX
#include <Windows.h>

class Timer {
private:
	LARGE_INTEGER nFreq_, nBefore_, nAfter_;
public:
	Timer() {}

	virtual ~Timer() {
	}

	void begin() {
		QueryPerformanceFrequency(&nFreq_);
		QueryPerformanceCounter(&nBefore_);
	}

	// ms
	float end() {
		QueryPerformanceCounter(&nAfter_);

		return ((nAfter_.QuadPart - nBefore_.QuadPart) * 1000.0f / nFreq_.QuadPart);
	}
};

#endif

} // namespace hstd

#endif // _TIMER_H_