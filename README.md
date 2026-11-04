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

For the moment the arborescence is like this : 


```text
simlab/
├── CMakeLists.txt        How to build the project; updated as the project grows
├── CONVENTIONS.md        All the conventions used in this project
├── src/
│   ├── core/             Elementary math tools needed throughout the project
│   ├── grid/             Grid data: Field2D, MacGrid2D, NodalGrid2D
│   ├── particles/        Particle storage (positions, velocities…): ParticleSet2D
│   ├── numerics/         Reusable numerical algorithms: advection, conjugate gradient, MAC Poisson…
│   ├── geometry/         Obstacles and collisions
│   ├── materials/        Physics specific to each material and its constitutive laws
│   ├── solvers/          Time loops producing the simulations: StableFluids, Flip, Mpm (all ISolver2D)
│   ├── scenes/           Initial conditions for the experiments
│   ├── diagnostics/      Validation quantities, exported to CSV
│   ├── render/           The only folder that uses OpenGL: all the graphics work
│   ├── gpu/              GPU implementation, if time allows
│   └── app/main.cpp      Everything in one place: the maestro!
├── shaders/              GLSL programs: *.vert, *.frag, *.comp
├── tests/                doctest unit tests: mini-checks that everything works
└── tools/                Python scripts to plot the results
```
