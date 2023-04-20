#pragma once
#include "raylib-cpp.hpp"

class Kernel {
protected:
	double supportRadius;
public:
	virtual double evaluate(double q) = 0;
	virtual raylib::Vector2 evaluateGradient(raylib::Vector2 x) = 0;
};