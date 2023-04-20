#pragma once
#include "raylib-cpp.hpp"
#include <vector>

class Particle {
public:
    Color color;

    int index = -1;

    double mass = 1;

    double radius = 2;

    raylib::Vector2 p;
    raylib::Vector2 v;

    raylib::Vector2 p0;
    raylib::Vector2 v0;

    double rho0;
    double rho;

    double aii0;
    double aii;

    raylib::Vector2 accViscosity;

    double psi;
    bool isBoundary;

    std::vector<Particle*> neighbours0;
    std::vector<Particle*> neighbours;

    double sourceTerm =0;
    raylib::Vector2 accPressure = raylib::Vector2::Zero();
    double Ap_i = 0;
    double pressure_DensitySolve = 0;
    double pressure_DivergenceSolve = 0;

    std::pair<int, int> cellKey = std::pair<int, int>(-1, -1);

    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    Particle(float x, float y, float vx, float vy, bool _boundary): p0(raylib::Vector2(x,y)), v0(raylib::Vector2(vx, vy)), isBoundary(_boundary) {
        reset();
        resetColor();
    }

    Particle(raylib::Vector2 _p, raylib::Vector2 _v, bool _boundary): p0(_p), v0(_v), v(raylib::Vector2::Zero()), isBoundary(_boundary) {
        reset();
        resetColor();
    }

    void resetColor() {
        if (isBoundary) {
            color = LIGHTGRAY;
        }
        else
        {
            color = BLUE;
        }
    }

    /**
     * Resets the position of this particle
     */
    void reset() {
        p = p0;
        v = v0;

        rho = rho0;
        aii = aii0;

        neighbours.clear();
        for (Particle * P: neighbours0) neighbours.push_back(P);
    }

    void setRest() {
        rho0 = rho;
        aii0 = aii;

        neighbours0.clear();
        for (Particle* P : neighbours) neighbours0.push_back(P);

        p0 = p;
        v0 = v;
    }
};