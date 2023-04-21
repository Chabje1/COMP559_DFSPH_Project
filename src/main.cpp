#include "raylib-cpp.hpp"
#include "FluidSPH.hpp"
#include "CubicSpline.hpp"
#include <string>

int screenWidth = 800;
int screenHeight = 450;


int main()
{
    raylib::Window window(screenWidth, screenHeight, "2D Divergence-Free SPH");
    raylib::Camera2D camera;

    camera.zoom = 1.0f;
    camera.offset = raylib::Vector2(screenWidth / 2, screenHeight / 2);
    camera.target = raylib::Vector2(screenWidth / 2, screenHeight / 2);
    camera.rotation = 0.0;

    FluidSPH<CubicSpline>* currentSim;
    std::vector<FluidSPH<CubicSpline>*> fluidSims;
    
    // =============== Simulation stuff ===============

    // Simulation with a rectangle of fluid falling on 2 lines
    FluidSPH<CubicSpline> fluid_sim1(0.001, &window, "Random Sim #1");

    fluid_sim1.createFluidRectangleFilled(15, 20, 30, 30);
    fluid_sim1.createBoundaryBox(fluid_sim1.particleRadius * fluid_sim1.scale, fluid_sim1.particleRadius * fluid_sim1.scale, 70, 70);
    fluid_sim1.createBoundaryLine(10, 200, 105, 265);
    fluid_sim1.createBoundaryLine(145, 265, 10, 345);

    fluidSims.push_back(&fluid_sim1);

    // Dam break
    FluidSPH<CubicSpline> fluid_sim2(0.001, &window, "Dam Break");

    fluid_sim2.createFluidRectangleFilled(10, 200, 30, 30);
    fluid_sim2.createBoundaryBox(fluid_sim2.particleRadius * fluid_sim2.scale, fluid_sim2.particleRadius * fluid_sim2.scale, 70, 70);

    fluidSims.push_back(&fluid_sim2);

    // Glass of Water
    FluidSPH<CubicSpline> fluid_sim3(0.001, &window, "Glass of Water");

    fluid_sim3.createFluidRectangleFilled(120, 150, 15, 15);
    fluid_sim3.createBoundaryBox(fluid_sim3.particleRadius * fluid_sim3.scale, fluid_sim3.particleRadius * fluid_sim3.scale, 70, 70);
    fluid_sim3.createBoundaryLine(125, 265, 145, 350);
    fluid_sim3.createBoundaryLine(165, 350, 185, 265);

    fluidSims.push_back(&fluid_sim3);

    // Dam break - Stress
    FluidSPH<CubicSpline> fluid_sim4(0.001, &window, "Dam Break - Stress");

    fluid_sim4.createFluidRectangleFilled(10, 50, 45, 79);
    fluid_sim4.createBoundaryBox(fluid_sim4.particleRadius * fluid_sim4.scale, fluid_sim4.particleRadius * fluid_sim4.scale, 100, 89);

    fluidSims.push_back(&fluid_sim4);

    // Dam break - Big Stress
    FluidSPH<CubicSpline> fluid_sim5(0.001, &window, "Dam Break - Big Stress");
    fluid_sim5.max_iter = 300;

    fluid_sim5.createFluidRectangleFilled(10, 50, 90, 79);
    fluid_sim5.createBoundaryBox(fluid_sim5.particleRadius * fluid_sim5.scale, fluid_sim5.particleRadius * fluid_sim5.scale, 159, 89);

    fluidSims.push_back(&fluid_sim5);

    //===============

    // Initialize the simulations
    for (FluidSPH<CubicSpline>* fsim : fluidSims) fsim->init();

    int currentSimInd = 0;
    currentSim = fluidSims[currentSimInd];

    // ================= Menu Stuff ===================

    int computationTime = 0;

    raylib::Vector2 helpMenuPosition(3. * screenWidth / 5., 0.);

    // The following are wrt to the help menu position
    float xoffset = 10;
    float yoffset = 10;
    raylib::Vector2 fpsPosition(xoffset, yoffset);

    raylib::Vector2 titlePosition(xoffset, fpsPosition.y + 20);
    int titleFontSize = 17;

    raylib::Vector2 dataPosition(xoffset, titlePosition.y + titleFontSize + 1);
    int dataFontSize = 14;
    
    std::vector<std::vector<std::pair<std::string, int*>>*> simDatas;
    std::vector<std::pair<std::string, int*>>* data;

    for (FluidSPH<CubicSpline>* sim : fluidSims) {

        data = new std::vector<std::pair<std::string, int*>>({
               {"Number of Particles :", &sim->numOfNonBoundary},
               {"Number of Iterations - Density :", &sim->numberOfIterations_Density},
               {"Number of Iterations - Divergence :", &sim->numberOfIterations_Divergence},
               {"Computation Time per Frame (ms): ", &computationTime}
        });

        simDatas.push_back(data);
    }

    raylib::Vector2 controlsPosition(xoffset, dataPosition.y + (dataFontSize + 1) * data->size());
    int controlsFontSize = 14;

    std::vector<std::string> controlsText = {
        "",
        "==== Controls ====",
        "Left-click: Select the Closest Particle",
        "Space: Pause/Unpause",
        "W/A/S/D: Move the camera",
        "Q/E: Zoom In/Out",
        "Z: Reset Camera",
        "R: Reset Simulation",
        "G: Draw Active Cells",
        "C: Color Particles wrt Velocities",
        "V: Draw Velocity Arrows",
        "F: Make Camera Follow the Selected Particle",
        "ESC: Quit",
        "O/P: Go to Next/Previous Simulation"
    };

    raylib::Vector2 helpMenuSize(screenWidth - helpMenuPosition.x, controlsPosition.y + (controlsFontSize + 1)*controlsText.size() + yoffset);

    bool isPaused = true;
    bool cameraFollowSelected = false;
    bool step = false;

    while (!window.ShouldClose())
    {
        // Selecting a particle to draw its neighbourhood
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            raylib::Vector2 mpos = camera.GetScreenToWorld(raylib::Mouse::GetPosition()) / currentSim->scale;
            Cell* cell = currentSim->spatialhash->GetCell(mpos);

            float closest = FLT_MAX;
            float dist;
            Particle* cl = NULL;
            for (Particle* p : cell->particles) {
                dist = (mpos - p->p).LengthSqr();
                if (dist < closest) {
                    closest = dist;
                    cl = p;
                }
            }

            if (currentSim->selected != NULL && cl != NULL && currentSim->selected->index == cl->index) {
                currentSim->selected = NULL;
            }
            else {
                currentSim->selected = cl;
            }
        }

        // Simulation Drawing and State Control
        if (IsKeyPressed(KEY_SPACE)) {
            isPaused = !isPaused;
        }
        else if (IsKeyPressed(KEY_R)) {
            currentSim->reset();
        }
        else if (IsKeyPressed(KEY_G)) {
            currentSim->drawCells = !currentSim->drawCells;
        }
        else if (IsKeyPressed(KEY_C)) {
            currentSim->colorVelocities = !currentSim->colorVelocities;
        }
        else if (IsKeyPressed(KEY_V)) {
            currentSim->drawVelocities = !currentSim->drawVelocities;
        }
        else if (IsKeyPressed(KEY_F)) {
            cameraFollowSelected = !cameraFollowSelected;
        }


        // Camera Position Control
        if (IsKeyDown(KEY_D)) {
            camera.target = Vector2Add(camera.target, Vector2(1, 0));
        }
        if (IsKeyDown(KEY_A)) {
            camera.target = Vector2Add(camera.target, Vector2(-1, 0));
        }
        if (IsKeyDown(KEY_W)) {
            camera.target = Vector2Add(camera.target, Vector2(0, -1));
        }
        if (IsKeyDown(KEY_S)) {
            camera.target = Vector2Add(camera.target, Vector2(0, 1));
        }

        if (currentSim->selected != NULL && cameraFollowSelected) {
            camera.target = currentSim->selected->p * currentSim->scale;
        }

        // Camera Zoom Control
        if (IsKeyPressed(KEY_Z)) {
            camera.zoom = 1.0f;
            camera.offset = raylib::Vector2(screenWidth / 2, screenHeight / 2);
            camera.target = raylib::Vector2(screenWidth / 2, screenHeight / 2);
            cameraFollowSelected = false;
        }
        if (IsKeyDown(KEY_Q)) {
            camera.zoom += 0.01f;
        }
        if (IsKeyDown(KEY_E)) {
            camera.zoom -= 0.01f;
        }

        if (camera.zoom < 0.1f) camera.zoom = 0.1f;

        if (IsKeyPressed(KEY_O)) {
            isPaused = true;
            if (++currentSimInd >= fluidSims.size()) currentSimInd = 0;
            currentSim = fluidSims[currentSimInd];
        }
        else if (IsKeyPressed(KEY_P)) {
            isPaused = true;
            if (--currentSimInd < 0) currentSimInd = fluidSims.size()-1;
            currentSim = fluidSims[currentSimInd];
        }

        if (IsKeyPressed(KEY_U)) {
            step = true;
        }

        BeginDrawing();
            ClearBackground(RAYWHITE);
            BeginMode2D(camera);
            currentSim->Draw();
            EndMode2D();
            
            // Draw Menu
            helpMenuPosition.DrawRectangle(helpMenuSize, Fade(SKYBLUE, 0.5f));
            DrawRectangleLines(helpMenuPosition.x, helpMenuPosition.y, helpMenuSize.x, helpMenuSize.y, BLUE);
            window.DrawFPS(helpMenuPosition.x + fpsPosition.x, helpMenuPosition.y + fpsPosition.y);

            raylib::DrawText(currentSim->name, helpMenuPosition.x + helpMenuSize.x/2 - dataFontSize / 2 *currentSim->name.length()/2, helpMenuPosition.y + titlePosition.y, titleFontSize, BLACK);

            int diter = 0;
            for (std::pair<std::string, int*> dataP : *simDatas[currentSimInd]) {
                raylib::DrawText(dataP.first, helpMenuPosition.x + dataPosition.x, helpMenuPosition.y + dataPosition.y + diter * (dataFontSize + 1), dataFontSize, BLACK);
                raylib::DrawText(std::to_string(*dataP.second), helpMenuPosition.x + dataPosition.x + dataFontSize/2 * dataP.first.length(), helpMenuPosition.y + dataPosition.y + diter * (dataFontSize + 1), dataFontSize, BLACK);
                diter++;
            }

            int cTiter = 0;
            for (std::string controltext : controlsText) {
                raylib::DrawText(controltext, helpMenuPosition.x + controlsPosition.x, helpMenuPosition.y + controlsPosition.y + cTiter * (controlsFontSize+1), controlsFontSize, BLACK);
                cTiter++;
            }

        EndDrawing();

        if (!isPaused || step) {
            computationTime = (int)std::floor(window.GetTime() * 1000);
            currentSim->Update();
            computationTime = (int)std::floor(window.GetTime() * 1000) - computationTime;
            step = false;
        }
    }

    return 0;
}