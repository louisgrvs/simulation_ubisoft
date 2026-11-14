# Conventions

## Numbers and Units
- `real = double` on the CPU (reference implementation); the GPU will use `float`.
- Arbitrary but consistent physical units: lengths in domain units, time in seconds.

## Coordinate System and Indexing
- Origin at the bottom-left, x pointing right, y pointing up (as in OpenGL).
- Cell `(i, j)` occupies `[i h, (i+1) h] x [j h, (j+1) h]`.
- Row-major storage, with `i` varying fastest: `k = i + nx * j`.
- Access through `f(i, j)`, never through `f[k]` outside `Field2D`.

## MAC Grid (`MacGrid2D`)

| Unknown | Dimensions | Sample Position `(i, j)` | Offset `(ox, oy)` |
|----------|------------|--------------------------|------------------|
| `p`, scalars | `nx x ny` | `((i+1/2) h, (j+1/2) h)` | `(1/2, 1/2)` |
| `u` | `(nx+1) x ny` | `(i h, (j+1/2) h)` | `(0, 1/2)` |
| `v` | `nx x (ny+1)` | `((i+1/2) h, j h)` | `(1/2, 0)` |
| Nodes (MPM) | `(nx+1) x (ny+1)` | `(i h, j h)` | `(0, 0)` |

- `u(i, j)` is the **left face** of cell `(i, j)`; `v(i, j)` is its **bottom face**.
- All sampling at arbitrary positions must go through `Field2D::sample` (bilinear interpolation, clamped at the boundaries).

## Walls (Days 4 and 5)
- Domain enclosed by 4 fixed walls with perfect slip (free-slip boundary conditions).
- Normal component: `u(0, j) = u(nx, j) = 0`, `v(i, 0) = v(i, ny) = 0` (`u·n = 0`).
  These faces are never treated as unknowns (in projection or diffusion).
- Tangential components and scalar fields (pressure, smoke): homogeneous Neumann boundary conditions; ghost cells outside the domain are simply omitted.

## Sign Conventions
- 2D vorticity: `\omega = \partial v/ \partial x - \partial u/\partial y`, positive in the counterclockwise direction.
- Pressure: solve `A q = -D w` with `A = -D G` (symmetric positive definite), `q = (dt/ρ) p`.
- Buoyancy: `f_y = \beta s`, pointing upward when `\beta > 0`.

## Rendering and Mouse Input
- `GL_R32F` texture: row 0 = `j = 0` = bottom of the domain; no vertical flipping.
- GLFW provides mouse coordinates with the origin at the TOP-left: invert y before converting to physical coordinates.

## Code
- Comments in French, identifiers in English (`camelCase` for functions, `PascalCase` for types, trailing `_` for private members).
- `src/`, excluding `render/` and `gpu/`: no dependencies on OpenGL.
- Every numerical mechanism must have a corresponding doctest test validating it.