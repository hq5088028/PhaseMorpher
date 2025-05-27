#pragma once
#include "../base/MACRO_DEF.h"

namespace pf {
	inline bool isTwoNumEquality(REAL a, REAL b) {
		if (std::fabs(a - b) < SYS_EPSILON)
			return true;
		else
			return false;
	}

	inline int REAL_to_int(REAL a) {
		if ((a - int(a)) > REAL(0.5))
			return int(a) + 1;
		else
			return int(a);
	}

	inline REAL interpolation_func(REAL phi) {
		return phi * phi * (3 - 2 * phi);
	}

	inline REAL dinterpolation_func(REAL phi) {
		return 6 * phi * (1 - phi);
	}

}