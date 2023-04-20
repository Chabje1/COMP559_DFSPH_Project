#pragma once
#include "Kernel.hpp"

class CubicSpline : public Kernel {
protected:
	double normCon;
    double normConGrad;
public:
	CubicSpline(double _supportRadius) {
		supportRadius = _supportRadius;
		normCon = 40.0 / (7.0 * (PI * _supportRadius * _supportRadius));
        normConGrad = 240.0 / (7.0 * (PI * _supportRadius * _supportRadius));
	}

	double evaluate(double q) {
		if (q > 1) return 0;
		else if (q <= 0.5) return normCon * (6.0 * pow(q, 3) - 6.0 * pow(q, 2) + 1.0);
		else return normCon * (2.0 * pow(1.0 - q, 3));
	}

	raylib::Vector2 evaluateGradient(raylib::Vector2 x) {
        double r = x.Length();
        double q = r / supportRadius;

        double v = 0;

        if (q <= 1) {
            if (r > 1e-6) {
                raylib::Vector2 grad = x / (r * supportRadius);
                if (q <= 0.5) {
                    return grad * normConGrad * q * (3.0 * q - 2.0);
                }
                else {
                    return grad * -normConGrad * pow(1.0 - q, 2);
                }
            }
        }
        return raylib::Vector2::Zero();
	}
};