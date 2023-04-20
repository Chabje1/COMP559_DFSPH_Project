#include "raylib-cpp.hpp"
#include "FluidSPH.hpp"
#include "CubicSpline.hpp"
#include <string>

int screenWidth = 800;
int screenHeight = 450;

int main()
{
    raylib::Window window(screenWidth, screenHeight, "Divergence-Free SPH");

    FluidSPH<CubicSpline> fluid(0.005, &window);

    fluid.createFluidRectangleFilled(15, 60, 30, 15);

    fluid.createBoundaryBox(fluid.particleRadius, fluid.particleRadius, 70, 70);
    fluid.createBoundaryLine(10, 200, 125, 275);
    fluid.createBoundaryLine(175, 300, 10, 400);

    fluid.init();

    while (!window.ShouldClose())
    {
        BeginDrawing();
            ClearBackground(RAYWHITE);
            fluid.Draw();
            window.DrawFPS();
            raylib::DrawText(std::to_string(fluid.numOfNonBoundary), 10, 30, 14, BLACK);
        EndDrawing();

        fluid.Update();
    }

    return 0;
}