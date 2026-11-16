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
- `Vec2` is a small home-made struct: it has `+`, `-`, `*` and `/` by a scalar, `+=`, `dot` and `norm`, but **no `-=`** (write `a = a - b`).
- Random numbers are drawn directly from `std::mt19937`, never through `std::uniform_real_distribution` (its algorithm differs between MSVC and GCC): same seed ⇒ same particles on every machine.

## Particles and Transfers (Day 7)
- Stratified seeding: each cell is split into `k x k` sub-cells (`k = 2` → 4 particles per cell), one point per sub-cell, with a jitter from 0 (regular lattice) to 1 (uniform in the sub-cell). A sub-cell is kept if its **centre** lies in the seeding region.
- Particle mass: `m_p = ρ h² / k²`, so the total mass is `ρ x seeded area`.
- P2G and G2P share a single stencil (`linearStencil`): the weights are exactly transposed.
  Exception to the `Field2D::sample` rule above: particle ↔ grid transfers and particle advection use this stencil (`interpolateMac`), so that G2P and the RK2 velocity use the same weights.
- **P2G ignores out-of-field nodes; G2P renormalises by the sum of the weights present** (constant extension of the tangential component).
- P2G does not enforce `u·n = 0`: that is the solver's job.

## Solids (Day 8)
- Signed distance `φ` (`SolidSdf2D`): **`φ > 0` in the free space** (water or air), **`φ < 0` inside a solid**, `φ = 0` on a wall. Normal `n = ∇φ / |∇φ|`, pointing from the solid towards the fluid.
- The tank walls are solids too: `φ = min(x, Lx − x, y, Ly − y, obstacles…)`. Nothing is stored on a grid: `φ` is computed on demand.
- On a MAC grid, solids are seen as a **staircase** of `Solid` cells (geometric error `O(h)`); particles see the true geometry through `φ`.

## Free Surface (Day 8)
- Cell labels: `Solid` if `φ(centre) < 0`, otherwise `Fluid` if the cell contains at least one particle, otherwise `Air`. Recomputed **at every step**.
- **Out of the domain, a cell is `Solid`**: the walls become Neumann conditions in the Laplacian with no special case. `enforceSolidBC` replaces `enforceWallBC` in the free-surface solver.
- **Isolated air pockets of at most 8 cells are filled** (→ `Fluid`). Otherwise, with `p = 0`, they implode under the hydrostatic pressure and the water volume melts away.
- Pressure unknowns live in `Fluid` cells only. **`p = 0` at the centre of `Air` cells** (Dirichlet); no term towards a `Solid` neighbour (Neumann). Row of a fluid cell:
  `(N q_ij − Σ_{Fluid neighbours} q_n) / h² = −div_ij`, with `N` = number of non-`Solid` neighbours. Non-fluid cells get an identity row (`q` stays 0).
- `A` is symmetric positive definite as soon as one `Air` cell exists. The mean of `q` is removed (`removeMean`) only when **all** cells are `Fluid` (pure Neumann).
- Consequence of `p = 0` at the air cell centre: the effective surface lies half a cell **above** the true one. Discrete hydrostatic reference: `p = ρ g (y_s − y)`, `y_s` = centre of the first air cell in the column.
- **A face is fluid** if one of its two cells is `Fluid` and none is `Solid`. The pressure gradient is subtracted on fluid faces only.
- A face touching a `Solid` cell has a zero normal component (`enforceSolidBC`); the tangential component is free (free slip).
- Velocity extrapolation: Jacobi-style layers (average of the valid neighbours, validated at the end of each layer), over `max(2, ⌈vmax dt / h⌉ + 2)` layers. Done twice per step: after P2G (faces with zero weight) and after the projection (from the fluid faces), each time followed by `enforceSolidBC`.
- Order of a step: labels → P2G → extrapolation → (`uOld`, `vOld`) → gravity and mouse → projection → extrapolation → G2P → advection.
- Gravity `g = (0, −9.81)` acts only in the `TankAtRest` scene, so the Day 7 scenes are unchanged.
- Advection: RK2 (midpoint) in the **grid** velocity, in `n` equal sub-steps, `n = ⌈vmax dt / (cfl h)⌉` with `cfl = 1` by default.
- After **every** sub-step, particles with `φ < margin` are moved back to `φ = margin` along the normal (up to 4 iterations for re-entrant corners), and their inward velocity is removed. `margin = 10⁻³ h`.
