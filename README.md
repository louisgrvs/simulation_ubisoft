# simulation_ubisoft
An introduction to simulation of different types of materials, using differents techniques.

This project is done by me, following Arnaud Schoenten's advice and the following articles : 
- https://graphicscomputing.fr/course/2024_2025/epita_ani3d/practice/content/05_fluids/a_stable_fluids/assets/stable_fluids.pdf
- https://www.cs.ubc.ca/~rbridson/docs/zhu-siggraph05-sandfluid.pdf
- https://disneyanimation.com/publications/the-affine-particle-in-cell-method/
- https://media.disneyanimation.com/uploads/production/publication_asset/94/asset/SSCTS13_2.pdf

The goal of this project is to do an introduction for  simulation, using different methods you can files in those articles. The main tool we will use will be OpenGL, as it does not ask for to many power from your computer and is pretty simple to learn :)

Everything will be in 2D, though a 3D implementation might be done depending of if I have the time or not.

The begining of this project will be consacrated to learning OpenGL following the link : https://learnopengl.com/Introduction.

We sometimes have sames files with ".h" or ".cpp" extension. The ".h" extension is made to declare objects (what exists and how do you call them) whereas the ".cpp" extension is made to implement all the methods that we need for a specific class. 

## Building the project

Requirements:
- CMake 3.16 or newer;
- a C++17 compiler (Visual Studio 2022 on Windows);
- Python 3 with `numpy` and `matplotlib`, for the plots;
- an internet connection for the first configuration only: CMake downloads GLFW and Dear ImGui automatically.

All the commands below are run from the folder that **contains** `simlab/` (the `build/` folder is created next to it):

```
cmake -S simlab -B build
cmake --build build --config Release
```

The first command configures the project (only needed once, or after moving the folder), the second one compiles everything. To compile a single program, add its name, for example:

```
cmake --build build --config Release --target advection_demo
```

All the executables end up in `build\Release\`.

## Running the programs

| Executable | Day | What it does | Window |
|---|---|---|---|
| `simlab_tests` | 2 to 7 | All the unit tests (verification code) | no |
| `advection_demo` | 3 | Zalesak's disk in a rotation: Euler vs RK2, several CFL numbers | no |
| `projection_demo` | 4 | Convergence of the conjugate gradient in the pressure projection | no |
| `smoke_demo` | 5 | Energy decay of a stationary vortex: numerical viscosity of Stable Fluids | no |
| `pic_transfer_demo` | 7 | Dissipation of the PIC transfers, and PIC vs Stable Fluids | no |
| `simlab` | 5 | Interactive 2D smoke (Stable Fluids) | yes |
| `simlab_particles` | 7 | Interactive 2D particles (PIC) | yes |

### Unit tests

At each step of our project, we use verification code (doctest). All the test files are compiled into a single program, `simlab_tests`. To run all the tests:

```
build\Release\simlab_tests.exe
```

To run only some of them:

```
build\Release\simlab_tests.exe -sf=*test_grid*           only one file
build\Release\simlab_tests.exe -tc=*Projection*          only the test cases whose name contains "Projection"
build\Release\simlab_tests.exe -ltc                      list all the test cases
build\Release\simlab_tests.exe -s -sf=*test_advection*   also show the checks that pass
```

### Console demos (CSV files + Python plots)

Each demo writes CSV files (and images for the advection) into an output folder. The matching Python script reads this folder and saves a figure in it.

```
build\Release\advection_demo.exe 128 sorties_advection
python simlab\tools\plot_advection.py sorties_advection

build\Release\projection_demo.exe sorties_projection
python simlab\tools\plot_projection.py sorties_projection

build\Release\smoke_demo.exe sorties_smoke
python simlab\tools\plot_smoke.py sorties_smoke

build\Release\pic_transfer_demo.exe sorties_pic
python simlab\tools\plot_pic.py sorties_pic sorties_smoke
```

- `advection_demo` takes the grid size first (here 128 × 128), then the output folder.
- The second folder given to `plot_pic.py` is optional: if it is the output of `smoke_demo`, the Stable Fluids curves are drawn next to the PIC ones.
- The figures are `advection.png`, `projection.png`, `smoke_energy.png` and `pic_transfers.png`, each in its output folder.

### Interactive applications

```
build\Release\simlab.exe
build\Release\simlab_particles.exe
```

They need an OpenGL 4.3 graphics driver. Their path to the `shaders/` folder is written at compile time, so they can be launched from any folder; if `simlab/` is moved, run the configuration command again.

Controls (both applications): `Space` pause, `N` or right arrow one step, `R` reset, `Esc` quit. Scenes, parameters and diagnostics are in the left panel.
- `simlab`: left click pushes the fluid and injects smoke, right click only injects smoke. The panel chooses which field is shown (smoke density, velocity norm, vorticity, pressure, divergence…).
- `simlab_particles`: left click pushes the fluid. The particles are drawn as small disks colored by their speed `|u_p|`; a grid field can be drawn behind them ("Fond" in the panel). In the "rotating disk" scene, the particles do not move: only the transfers between particles and grid are applied, to measure how much PIC dissipates. Keep the color scale fixed to see their speed fade.

### From Visual Studio

Right click on a program, "Set as Startup Project", then `Ctrl+F5`. The arguments of a console demo go in the project properties, "Debugging" -> "Command Arguments" (for example `128 sorties_advection`).

## How the display works (OpenGL)

The simulation and the display are fully separated:
- the simulation (`simcore`: `core/`, `grid/`, `numerics/`, `particles/`, `scenes/`, `diagnostics/`, `solvers/`) runs on the CPU, in `double`, and never calls OpenGL;
- the display (`render/`, `shaders/`, and the two `main` files in `app/`) only **reads** the simulation's data and turns it into pixels on the GPU, in `float`.

So a bug in a picture is either in the data (check it with the tests and the console demos) or in the display, never in both at once.

### What happens at each frame

The two interactive applications (`app/main.cpp` and `app/particles_main.cpp`) follow the same loop:

1. **Window and OpenGL context.** At start-up, GLFW opens the window and creates an OpenGL 4.3 Core context; glad then loads the OpenGL functions (`external/glad/`); ImGui creates the panel.
2. **Events.** Keyboard keys, and the mouse position converted from pixels to physical coordinates (meters) to drive the brush.
3. **Simulation step.** `solver.step(dt)`, unless the simulation is paused.
4. **Rendering.** The data are sent to the GPU and drawn (see below), in a rectangle on the right of the panel that keeps the proportions of the domain.
5. **Panel.** ImGui draws the panel on top, then the image is shown (`glfwSwapBuffers`).

### Fields: `FieldRenderer` + `field.vert` / `field.frag`

```
Field2D (CPU, double) --FieldRenderer::upload--> texture GL_R32F (GPU, float) --field.vert + field.frag--> colored rectangle
```

- `FieldRenderer::upload` copies a `Field2D` (smoke, speed, vorticity, pressure…) into a one-channel float texture: **one texel per cell**.
- `field.vert` (vertex shader: runs once per vertex, decides *where* things are drawn) draws a rectangle covering the viewport. Its 4 corners are computed from `gl_VertexID`, so no vertex buffer is needed.
- `field.frag` (fragment shader: runs once per pixel, decides *its color*) reads the texture at this pixel and turns the value into a color: the "inferno" palette (black → yellow) for positive fields, a diverging palette (blue < 0, orange > 0) for signed fields such as vorticity or pressure. NaN or infinite values are drawn in magenta, so a numerical explosion is visible at once.
- Interpolation between cells is done by the GPU itself (texture filtering): smooth (bilinear) or blocky (one color per cell), switchable in the panel.

### Particles: `ParticleRenderer` + `particles.vert` / `particles.frag`

```
ParticleSet2D (CPU, double) --ParticleRenderer::upload--> vertex buffer [x, y, |u|] (GPU, float) --particles.vert + particles.frag--> disks
```

- `ParticleRenderer::upload` writes, for each particle, its position and its speed into a vertex buffer. The buffer is allocated once, then only updated at each frame (`glBufferSubData`).
- Each particle is drawn as a point (`GL_POINTS`). `particles.vert` converts its physical position into screen coordinates (`2 x / L - 1`) and sets its size in pixels (`gl_PointSize`).
- `particles.frag` receives a square of pixels for each point and keeps only the disk inside it (`gl_PointCoord`); the color comes from the same "inferno" palette as the fields.
- The field (if any) is drawn first, the particles on top, in the same rectangle: the two always match.

### Modifying the display

- The shaders are text files read **when the application starts**, they are not compiled by CMake: after editing a `.vert` or `.frag`, just restart the application, no need to rebuild.
- If a shader does not compile, the application stops and prints the GLSL error message with the line number.
- OpenGL errors are printed in the console (debug callback enabled in the `main` files).
- To understand these files, the useful LearnOpenGL chapters are *Hello Triangle*, *Shaders*, *Textures* and *Advanced Data*.

## Project tree

For the moment the arborescence is like this : 

```text
simlab/
├── CMakeLists.txt        How to build the project; updated as the project grows
├── CONVENTIONS.md        All the conventions used in this project
├── external/             Third-party code shipped with the project
│   └── glad/             OpenGL 4.3 Core function loader (generated, never edited by hand)
│       ├── include/
│       │   ├── glad/
│       │   │   └── glad.h           Declarations of all the OpenGL 4.3 functions and constants
│       │   └── KHR/
│       │       └── khrplatform.h    Platform types needed by glad.h
│       └── src/
│           └── glad.c               Loads the address of every OpenGL function from the graphics driver at start-up
├── src/
│   ├── core/             Elementary math tools needed throughout the project
│   │   ├── Types.h                  Obsolete since day 5 (the type `real` is now in Vec2.h): no longer compiled, can be deleted
│   │   ├── Vec2.h                   Definition of some useful type in this project, and operation on those types (for example 2D vectors)
│   │   └── Kernels.h                #TODO (jour 7)
│   ├── grid/             Grid data: Field2D, MacGrid2D, NodalGrid2D
│   │   ├── Field2D.h                A class of a basic grid on which we are going to interpolate everything, a MACgrid basically consist in 3 field2D grid
│   │   ├── Field2D.cpp              Methode for this class
│   │   ├── MacGrid2D.h              A class that will consitute our MacGrid class, constructed as 3 field2D.h elements
│   │   └── MacGrid2D.cpp            Methods for this class
│   ├── particles/        Particle storage (positions, velocities…): ParticleSet2D
│   │   └── ParticleSet2D.h          #TODO (jour 7)
│   ├── numerics/         Reusable numerical algorithms: advection, conjugate gradient, MAC Poisson…
│   │   ├── Advection.h              Advection method, as described in STAM articles. We use the caracteristic method to get our approximation of NVS solution
│   │   ├── Advection.cpp            Method for the advection.h class.
│   │   ├── MacOperators.h           Class of diferential operator : laplacian, divergence, gradient ...
│   │   ├── MacOperators.cpp         Method to compute those operators
│   │   ├── ConjugateGradient.h      A conjugate gradient algorithme useful to solve equation of type Ax = y. Useful for the projection algorithm
│   │   ├── ConjugateGradient.cpp    Method to compute the algorithme
│   │   ├── Projection.h             Projection's method following STAM (how to obtain w_4)
│   │   ├── Projection.cpp           Function to compute this method
│   │   ├── Diffusion.h              Diffusion's method following STAM (how to obtain w_3)
│   │   ├── Diffusion.cpp            Function to compute this method
│   │   ├── Forces.h                 Force's method follow STAM (how to obtain w_1)
│   │   ├── Forces.cpp               Function to compute this method
│   │   ├── MacTransfer.h            #TODO (jour 7)
│   │   └── MacTransfer.cpp          #TODO (jour 7)
│   ├── geometry/         Obstacles and collisions
│   ├── materials/        Physics specific to each material and its constitutive laws
│   ├── solvers/          Time loops producing the simulations: StableFluids, Flip, Mpm (all ISolver2D)
│   │   ├── ISolver2D.h              A class that will simulate our animation, using time dependency, 
│   │   ├── TestPatternSolver.h      Obsolete since day 5 (test pattern of day 2): no longer compiled, can be deleted
│   │   ├── TestPatternSolver.cpp    Obsolete since day 5 (test pattern of day 2): no longer compiled, can be deleted
│   │   ├── StableFluidsSolver.h     A solver as given in description in the STAM articles: one time step chains forces (w_1), advection (w_2), diffusion (w_3) and projection (w_4)
│   │   ├── StableFluidsSolver.cpp   Methods for this solver
│   │   ├── FlipSolver.h             A solver as given in description in the zhu-brightson article, it is a mixt between FLIP and PIC methods (PIC only for now)   #TODO (jour 7)
│   │   └── FlipSolver.cpp           Methods for this solver   #TODO (jour 7)
│   ├── scenes/           Initial conditions for the experiments
│   │   ├── RotationScene.h          Some definition of useful function we use to check that our code works well. We use a basic rotation field and the Zalesak's disk
│   │   ├── SmokeScenes.h            A stationnary Euler solution that we will use to simulate our smoke : it will be the initial configuration
│   │   ├── Seeding.h                #TODO (jour 7)
│   │   └── Seeding.cpp              #TODO (jour 7)
│   ├── diagnostics/      Validation quantities, exported to CSV
│   │   ├── FluidDiagnostics.h       Quantities that we will use to check that our simulation deos what we need
│   │   ├── FluidDiagnostics.cpp     Actual implementation of those quantities
│   │   ├── ParticleDiagnostics.h    #TODO (jour 7)
│   │   └── ParticleDiagnostics.cpp  #TODO (jour 7)
│   ├── render/           The only folder that uses OpenGL: all the graphics work
│   │   ├── FieldRenderer.h          Draws a Field2D as a colored image: one texel per cell, palette chosen in field.frag
│   │   ├── FieldRenderer.cpp        Compiles field.vert/field.frag, uploads the field into a float texture, draws it in a rectangle
│   │   ├── ParticleRenderer.h       Draws the particles as colored disks
│   │   └── ParticleRenderer.cpp     Compiles particles.vert/particles.frag, fills the vertex buffer (position, speed), draws GL_POINTS
│   ├── gpu/              GPU implementation, if time allows
│   └── app/              Executables
│       ├── main.cpp                 Everything in one place: the maestro! Interactive smoke application (window, mouse, solver, display, panel)
│       ├── advection_demo.cpp       #TODO (jour 3)
│       ├── projection_demo.cpp      #TODO (jour 4)
│       ├── smoke_demo.cpp           #TODO (jour 5)
│       ├── particles_main.cpp       Interactive particles application (PIC): same loop as main.cpp, particles drawn on top of a grid field
│       └── pic_transfer_demo.cpp    #TODO (jour 7)
├── shaders/              GLSL programs run by the GPU: *.vert, *.frag, *.comp (read at start-up, no rebuild needed)
│   ├── field.vert                   Vertex shader: a rectangle covering the viewport, corners computed from gl_VertexID
│   ├── field.frag                   Fragment shader: reads the field texture and maps each value to a color (inferno or diverging palette, NaN in magenta)
│   ├── particles.vert               Vertex shader: physical position -> screen coordinates, point size in pixels
│   └── particles.frag               Fragment shader: keeps a disk in each point, colored by the particle speed
├── tests/                doctest unit tests: mini-checks that everything works
│   ├── doctest.h                    #TODO (jour 2)
│   ├── test_main.cpp                #TODO (jour 2)
│   ├── test_field2d.cpp             Obsolete since day 5 (replaced by test_grid.cpp): no longer compiled, can be deleted
│   ├── test_grid.cpp                
│   ├── test_advection.cpp           #TODO (jour 3)
│   ├── test_projection.cpp          #TODO (jour 4)
│   ├── test_smoke.cpp               #TODO (jour 5)
│   └── test_pic.cpp                 #TODO (jour 7)
└── tools/                Python scripts to plot the results
    ├── plot_advection.py            A python code to show how our advection method (from advection.cpp) works, using Zelestak's disk. You can check the results in sorties_advection
    ├── plot_projection.py           #TODO (jour 4)
    ├── plot_smoke.py                #TODO (jour 5)
    └── plot_pic.py                  #TODO (jour 7)
```
