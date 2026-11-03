# simulation_ubisoft
An introduction to real-time simulation of different types of materials

This project is done by me, following Arnaud Schoenten's advice and the following articles : 
- https://graphicscomputing.fr/course/2024_2025/epita_ani3d/practice/content/05_fluids/a_stable_fluids/assets/stable_fluids.pdf
- https://www.cs.ubc.ca/~rbridson/docs/zhu-siggraph05-sandfluid.pdf
- https://disneyanimation.com/publications/the-affine-particle-in-cell-method/
- https://media.disneyanimation.com/uploads/production/publication_asset/94/asset/SSCTS13_2.pdf

The goal of this project is to do an introduction for a real-time simulation, using different methods you can files in those articles. The main tool we will use will be OpenGL, as it does not ask for to many power from your computer :)

Everything will be in 2D, though a 3D implementation might be done depending of if I have the time or not.

The begining of this project will be consacrated to learning OpenGL following the link : https://learnopengl.com/Introduction.

For the moment the arborescence is like this : 
simlab/
├── CMakeLists.txt            It will describe how to compile the project, I will ugprade it following how I advance.
├── CONVENTIONS.md            All the conventions I will eb using in this project.
├── src/
│   ├── core/                 Every elementary mathematics tools that I would need for the project.
│   ├── grid/                 Grid data : Field2D, MacGrid2D, NodalGrid2D
│   ├── particles/            particle storage : positions, speeds etc ... ParticleSet2D
│   ├── numerics/             numercis algrotihm re-usable : advection, conjuguated gradient , Poisson MAC...
│   ├── geometry/             obstacles and collisions : 
|   ├── materials/            physics proper to each materials describe and the laws that implemented it
│   ├── solvers/              Time loop to obtain our simulation : StableFluids, Flip, Mpm (all ISolver2D)
│   ├── scenes/               All of our initial conditions for our experiments
│   ├── diagnostics/          All of our validations things
│   ├── render/               The only folder that uses OpenGL : it is where our graphics works.
│   ├── gpu/                  If we have time : the GPU implementation
│   └── app/main.cpp          Everything in one : the maestro !
├── shaders/                  GLSL's programs *.vert, *.frag, *.comp
├── tests/                    doctest : mini-tests to check we're doing allright !
└── tools/                    scripts Python to draw our results.
