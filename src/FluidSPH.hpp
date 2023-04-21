#pragma once
#include "Particle.hpp"
#include "Kernel.hpp"
#include "SpatialHash.hpp"
#include <vector>

/**
* A class that implements the fluid simulation algorithm DFSPH.
* 
* - The type T must be a class who has Kernel as its base class.
* An example of this is the CubicSpline class.
*/
template <class T> requires std::is_base_of<Kernel, T>::value class FluidSPH {
public:
    // Simulation Settings
    bool useGravity = true;
    raylib::Vector2 gravity = raylib::Vector2(0, 9.81);
    bool drawCells = false;
    bool colorVelocities = false;
    bool drawVelocities = false;
    float velocityScale = 0.1;

    std::string name;

    // Fluid Properties
	double viscousDamping = 0.05;
    double material_density = 1000;

    // Particle Properties
    double particleRadius = 0.025;
    double supportRadius = particleRadius * 4;
	double neighbourhoodDistanceSquared = pow(supportRadius,2);
    double pMass = 4*pow(particleRadius, 2) * material_density;

    float scale = 100;

    // Solver Settings
    double omega = 0.5;
    int max_iter = 100;
    double maxError_Density = 0.01;
    double maxError_Divergence = 0.1;

    // =================

    int numberOfIterations_Density = 0;
    int numberOfIterations_Divergence = 0;

    double dt;
    int next_part_id = 0;
    int numOfNonBoundary = 0;

    raylib::Window* parent;

    Kernel* kernel = new T(supportRadius);

    std::vector<Particle*> particles;
    std::vector<Particle*> nonBoundaryParticles;

    SpatialHash* spatialhash;

    // ================

    Particle* selected = NULL;


    // Constructor
    FluidSPH(double timestep, raylib::Window* _parent, std::string _name): name(_name) {
        dt = timestep;
        parent = _parent;
        spatialhash = new SpatialHash(supportRadius, 140*2* supportRadius, 140 * 2 * supportRadius);
	}

	/**
	* Resets the fluid to its initial state
	*/
	void reset() {
        for (Particle* p : nonBoundaryParticles) p->reset();
        for (Particle* p : nonBoundaryParticles) spatialhash->AddParticleToCell(p);
	}

    // ========================================================================
    //                            Fluid Creating Stuff 
    // ========================================================================

    /**
    Creates a new particle in the fluid at position (x,y) with initial velocity (vx,vy).
    isBoundary sets whether the particle is a boundary particle. (ie. static)
    */
    Particle* createParticle(float x, float y, float vx, float vy, bool isBoundary) {
        Particle* p = new Particle(x, y, vx, vy, isBoundary);

        p->radius = particleRadius;
        p->mass = pMass;
        particles.push_back(p);
        p->index = next_part_id++;

        if (!isBoundary) {
            nonBoundaryParticles.push_back(p);
            numOfNonBoundary++;
        }

        spatialhash->AddParticleToCell(p);
        return p;
    }

    /**
    Creates a new particle in the fluid at position _p with initial velocity _v.
    isBoundary sets whether the particle is a boundary particle. (ie. static)
    */
    Particle* createParticleV(raylib::Vector2 _p, raylib::Vector2 _v, bool isBoundary) {
        Particle* p = new Particle(_p, _v, isBoundary);
        p->radius = particleRadius;
        p->mass = pMass;
        particles.push_back(p);
        p->index = next_part_id++;
        if (!isBoundary) {
            nonBoundaryParticles.push_back(p);
            numOfNonBoundary++;
        }

        spatialhash->AddParticleToCell(p);
        return p;
    }

    /**
    Creates a rectangle of fluid with its top-right corner at (x1,y1) and its bottom-left corner at (x2,y2). 
     xSpacing & ySpacing control the spacing between the particles' centers.
     This means that a spacing of 2 * particleRadius has them touching.
    */
    void createFluidRectangle(float _x1, float _y1, float _x2, float _y2, float _xSpacing, float _ySpacing) {
        float x1 = _x1 / scale; float y1 = _y1 / scale; float x2 = _x2 / scale; float y2 = _y2 / scale; float xSpacing = _xSpacing / scale; float ySpacing = _ySpacing / scale;
        for (float i = x1; i <= x2; i += xSpacing) {
            for (float j = y1; j <= y2; j += ySpacing) {
                createParticle(i, j, 0, 0, false);
            }
        }
    }

    /**
    Creates a completely filled rectangle of fluid which is numX particles wide and numY particles high
    */
    void createFluidRectangleFilled(float _x, float _y, int numX, int numY) {
        float x = _x / scale; float y = _y / scale;
        for (float i = 0; i < numX; i++) {
            for (float j = 0; j < numY; j++) {
                createParticle(i * 2*particleRadius + x , 2*j * particleRadius + y, 0, 0, false);
            }
        }
    }

    /**
    Creates a Box of Boundary particles with its top-right corner at (x1,y1) and is numX wide and numY high
    */
    void createBoundaryBox(float _x, float _y, int numX, int numY) {
        float x = _x / scale; float y = _y / scale;
        for (float i = 1; i < numX; i++) {
            createParticle(x + i * 2 * particleRadius, y, 0, 0, true);
            createParticle(x + i * 2 * particleRadius, y + numY * 2 * particleRadius, 0, 0, true);
        }

        for (float i = 1; i < numY; i++) {
            createParticle(x,y + i * 2 * particleRadius, 0, 0, true);
            createParticle(x + numX * 2 * particleRadius, y + i * 2 * particleRadius, 0, 0, true);
        }

        createParticle(x, y, 0, 0, true);
        createParticle(x, y + numY * 2 * particleRadius, 0, 0, true);
        createParticle(x + numX * 2 * particleRadius, y, 0, 0, true);
        createParticle(x + numX * 2 * particleRadius, y + numY * 2 * particleRadius, 0, 0, true);
    }

    /**
    Creates a Line of Boundary particles with from (x1,y1) to (x2,y2).
    */
    void createBoundaryLine(float _x1, float _y1, float _x2, float _y2) {
        float x1 = _x1 / scale; float y1 = _y1 / scale; float x2 = _x2 / scale; float y2 = _y2 / scale;
        raylib::Vector2 p = raylib::Vector2(x1, y1);
        raylib::Vector2 dir = raylib::Vector2(x2, y2) - p;

        double dirLen = dir.Length();
        raylib::Vector2 dp = dir/dirLen * 2 * particleRadius;

        for (int i = 0; i * 2 * particleRadius < dirLen; i++) {
            createParticleV(p, raylib::Vector2::Zero(), true);
            p += dp;
        }
    }

    // ========================================================================
    //                            Fluid Update Stuff 
    // ========================================================================

    /**
    Finds the neighbours of Particle p by checking its cell and those that neighbour it.
    Its cell is assigned to it by the SpatialHash class.
    */
    void findNeighbours(Particle* p) {
        p->neighbours.clear();

        Cell* particleCell = spatialhash->GetCell(p);

        for (Particle* o : particleCell->particles) {
            if (p->index != o->index && (p->p - o->p).LengthSqr() <= neighbourhoodDistanceSquared) {
                p->neighbours.push_back(o);
            }
        }

        for (Cell* cell : particleCell->neighbours) {
            for (Particle* o : cell->particles) {
                if (p->index != o->index && (p->p - o->p).LengthSqr() <= neighbourhoodDistanceSquared) {
                    p->neighbours.push_back(o);
                }
            }
        }
    }
    
    /**
    Initializes the simulation by computing all the starting values for 
    the neighbourhoods, the densities, the A_ii's and the psi's.
    */
    void init() {
        for (Particle* pi : particles) {
            findNeighbours(pi);

            if (!pi->isBoundary) {
                UpdateNonBoundaryParticleProperties(pi);

                pi->setRest();
            }
            else ComputePsi(pi);
        }
    }

    /**
    Computes the neighbourhoods, the densities and the A_ii's for non-boundary Particle pi.
    */
    void UpdateNonBoundaryParticleProperties(Particle* pi) {
        findNeighbours(pi);

        pi->rho = pi->mass * kernel->evaluate(0);

        for (Particle* pj : pi->neighbours)
        {
            if (pj->isBoundary) {
                pi->rho += pj->psi * kernel->evaluate((pi->p - pj->p).Length() / supportRadius);
            }
            else {
                pi->rho += pj->mass * kernel->evaluate((pi->p - pj->p).Length() / supportRadius);
            }
        }


        raylib::Vector2 a = raylib::Vector2::Zero();
        double b = 0;
        raylib::Vector2 buf;

        for (Particle* pj : pi->neighbours)
        {
            if (!pj->isBoundary) {
                buf = kernel->evaluateGradient(pi->p - pj->p) * pj->mass;
            }
            else {
                buf = kernel->evaluateGradient(pi->p - pj->p) * pj->psi;
            }
            a += buf;
            b += pow(buf.Length(), 2);
        }
        pi->aii = -dt / pow(pi->rho, 2) * (pow(a.Length(), 2) + b);
    }

    /**
    Computes the pseudo-mass value for a given boundary particle pi.
    */
    void ComputePsi(Particle* pi) {
        pi->psi = 0;
        for (Particle* pj : pi->neighbours)
        {
            if (pj->isBoundary) pi->psi += kernel->evaluate((pi->p - pj->p).Length() / supportRadius);
        }
        pi->psi = material_density / pi->psi;
    }

    /**
    Steps the simulation by dt.
    */
    void Update() {
        // VISCOSITY
        for (Particle* pi : nonBoundaryParticles) {
            pi->accViscosity = raylib::Vector2::Zero();

            for (Particle* pj : pi->neighbours) {
                if (!pj->isBoundary) {
                    pi->accViscosity -= (pi->v - pj->v) * pj->mass / pj->rho * kernel->evaluate((pi->p - pj->p).Length() / supportRadius);
                }
            }
            pi->accViscosity *= viscousDamping / dt;
        }

        // INTEGRATE VELOCITIES WITH NON-PRESSURE FORCES
        for (Particle* pi : nonBoundaryParticles) {
            pi->v += (pi->accViscosity + ((useGravity)?gravity:raylib::Vector2::Zero())) * dt;
        }

        // COMPUTE KAPPA TERMS
        double pi_star = 0;
        for (Particle* pi : nonBoundaryParticles) {
            pi_star = pi->rho;
            for (Particle* pj : pi->neighbours) {
                if (pj->isBoundary) {
                    pi_star += dt * pj->psi * (pi->v).DotProduct(kernel->evaluateGradient(pi->p - pj->p));
                }
                else {
                    pi_star += dt * pj->mass * (pi->v - pj->v).DotProduct(kernel->evaluateGradient(pi->p - pi->p));
                }
            }

            pi->sourceTerm = (material_density - pi_star) / dt;
        }

        double avg_density_error = 10000000;
        double density_error;
        numberOfIterations_Density = 0;

        for (Particle* pi : nonBoundaryParticles) {
            pi->pressure_DensitySolve = 0;
            pi->pressure_DivergenceSolve = 0;
        }

        while ((avg_density_error > maxError_Density * 0.01 * material_density && max_iter > numberOfIterations_Density) || numberOfIterations_Density < 2) {
            // Compute Pressure Acceleration
            for (Particle* pi : nonBoundaryParticles) {
                pi->accPressure = raylib::Vector2::Zero();

                for (Particle* pj : pi->neighbours) {
                    if (pj->isBoundary) {
                        pi->accPressure -= kernel->evaluateGradient(pi->p - pj->p) *pj->psi* pi->pressure_DensitySolve / pow(pi->rho, 2);
                    }
                    else {
                        pi->accPressure -= kernel->evaluateGradient(pi->p - pj->p) *pj->mass* (pi->pressure_DensitySolve / pow(pi->rho, 2) + pj->pressure_DensitySolve / pow(pj->rho, 2));
                    }
                }
            }

            // Compute (Ap)_i
            for (Particle* pi : nonBoundaryParticles) {
                pi->Ap_i = 0;

                for (Particle* pj : pi->neighbours) {
                    if (pj->isBoundary) {
                        pi->Ap_i += pj->psi * pi->accPressure.DotProduct(kernel->evaluateGradient(pi->p - pj->p));
                    }
                    else {
                        pi->Ap_i += pj->mass * (pi->accPressure - pj->accPressure).DotProduct(kernel->evaluateGradient(pi->p - pj->p));
                    }
                }

                pi->Ap_i *= dt;
            }

            density_error = 0;
            // Update Pressure Values
            for (Particle* pi : nonBoundaryParticles) {
                if (std::abs(pi->aii) > 1e-6) {
                    pi->pressure_DensitySolve += omega / pi->aii * (pi->sourceTerm - pi->Ap_i);
                }
                else pi->pressure_DensitySolve = 0;
                pi->pressure_DensitySolve = std::fmax(pi->pressure_DensitySolve, 0);
                density_error -= std::fmin(pi->sourceTerm - pi->Ap_i, 0) * dt;
            }
            avg_density_error = density_error / numOfNonBoundary;
            numberOfIterations_Density++;
        }

        // INTEGRATE VELOCITIES & POSITIONS WITH PRESSURE FORCE
        for (Particle* pi : nonBoundaryParticles) {
            pi->v += pi->accPressure * dt;
            pi->p += pi->v * dt;
            spatialhash->AddParticleToCell(pi);
        }

        // Recompute neighbourhoods, densities and aii
        for (Particle* pi : nonBoundaryParticles) {
            UpdateNonBoundaryParticleProperties(pi);
        }

        // COMPUTE DIVERGENCE SOURCE TERM
        for (Particle* pi : nonBoundaryParticles) {
            pi->sourceTerm = 0;
            for (Particle* pj : pi->neighbours) {
                if (pj->isBoundary) {
                    pi->sourceTerm -= pj->psi * (pi->v).DotProduct(kernel->evaluateGradient(pi->p - pj->p));
                }
                else {
                    pi->sourceTerm -= pj->mass * (pi->v - pj->v).DotProduct(kernel->evaluateGradient(pi->p - pj->p));
                }
            }
        }

        // Do the Divergence Solve
        avg_density_error = 10000000;
        numberOfIterations_Divergence = 0;

        while ((avg_density_error > maxError_Divergence * 0.01 * material_density && max_iter > numberOfIterations_Divergence) || numberOfIterations_Divergence < 2) {
            // Compute Pressure Acceleration
            for (Particle* pi : nonBoundaryParticles) {
                pi->accPressure = raylib::Vector2::Zero();

                for (Particle* pj : pi->neighbours) {
                    if (pj->isBoundary) {
                        pi->accPressure -= kernel->evaluateGradient(pi->p - pj->p) * pj->psi * pi->pressure_DivergenceSolve / pow(pi->rho, 2);
                    }
                    else {
                        pi->accPressure -= kernel->evaluateGradient(pi->p - pj->p) * pj->mass * (pi->pressure_DivergenceSolve / pow(pi->rho, 2) + pj->pressure_DivergenceSolve / pow(pj->rho, 2));
                    }
                }
            }

            // Compute (Ap)_i
            for (Particle* pi : nonBoundaryParticles) {
                pi->Ap_i = 0;

                for (Particle* pj : pi->neighbours) {
                    if (pj->isBoundary) {
                        pi->Ap_i += pj->psi * pi->accPressure.DotProduct(kernel->evaluateGradient(pi->p - pj->p));
                    }
                    else {
                        pi->Ap_i += pj->mass * (pi->accPressure - pj->accPressure).DotProduct(kernel->evaluateGradient(pi->p - pj->p));
                    }
                }

                pi->Ap_i *= dt;
            }


            density_error = 0;
            // Update Pressure Values
            for (Particle* pi : nonBoundaryParticles) {
                if (std::abs(pi->aii) > 1e-6) {
                    pi->pressure_DivergenceSolve += omega / pi->aii * (pi->sourceTerm - pi->Ap_i);
                }
                else pi->pressure_DivergenceSolve = 0;
                pi->pressure_DivergenceSolve = std::fmax(pi->pressure_DivergenceSolve, 0);
                density_error -= std::fmin(pi->sourceTerm - pi->Ap_i, 0) * dt;
            }
            avg_density_error = density_error / numOfNonBoundary;
            numberOfIterations_Divergence++;
        }

        // INTEGRATE VELOCITIES WITH PRESSURE FORCE
        for (Particle* pi : nonBoundaryParticles) {
            pi->v += pi->accPressure * dt;
        }

    }

	void Draw() {
        if (drawCells) {
            for (Particle* p : particles)
            {
                DrawSpatialHashCell(p->cellKey);
            }
        }

        for (Particle* p : particles) {
            if (colorVelocities && !p->isBoundary) {
                double vl = p->v.LengthSqr()/scale;
                double sig = vl / (1 + std::abs(vl)) * 255;
                p->color = raylib::Color(sig, sig, 255-sig);
            }
            else p->resetColor();

            if (drawVelocities && !p->isBoundary) {
                (p->p * scale).DrawLine((p->p + p->v * velocityScale)*scale);
            }
        }

        if (selected != NULL) {
            selected->color = RED;
            for (Particle* p : selected->neighbours) {
                p->color = GREEN;
            }
        }

		for (Particle* p : particles)
		{
			(p->p * scale).DrawCircle((p->radius * scale), p->color);
		}
	}

    void DrawSpatialHashCell(std::pair<int, int> k) {
        Cell* c = spatialhash->GetCell(k);
        c->Draw(scale);
    }

    void DrawSpatialHashCellandNeighbours(int i, int j) {
        Cell* c = spatialhash->GetCell(std::pair<int, int>(i, j));
        c->Draw(scale);
        c->DrawNeighbours(scale);
    }
};