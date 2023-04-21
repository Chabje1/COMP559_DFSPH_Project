#pragma once
#include "Kernel.hpp"

class CubicSpline : public Kernel {
protected:
	double sigma;
public:
	CubicSpline(double _supportRadius) {
		supportRadius = _supportRadius;
		sigma = 40.0 / (7.0 * (PI * _supportRadius * _supportRadius));
	}

	double evaluate(double q) {
		// Cubic spline equations if we scale the q by 2
		// since we want it to be only non-zero for 
		// 0 <= q <= 1. In the original, Cubic Spline,
		// the function is non-zero for 0 <= q <= 2.
		if (q > 1) return 0;
		else if (q <= 0.5) return sigma * (6.0 * pow(q, 3) - 6.0 * pow(q, 2) + 1.0);
		else return sigma * (2.0 * pow(1.0 - q, 3));
	}

	raylib::Vector2 evaluateGradient(raylib::Vector2 x) {
        double r = x.Length();
		if (r < 1e-6) return raylib::Vector2::Zero();
        double q = r / supportRadius;

		// Grad W(q(x)) = x/(h*||x||) * dWdq 
        if (q > 1) return raylib::Vector2::Zero();
		else if (q <= 0.5) return x / (r * supportRadius) * 6 * sigma * (3.0 * pow(q,2) - 2.0 * q);
		else return x / (r * supportRadius) * -6 * sigma * pow(1.0 - q, 2);
        
	}
};