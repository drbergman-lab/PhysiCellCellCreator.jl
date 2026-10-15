# Design brief: expanding the patch shapes (ellipses now, 3D next)

Date: 2026-10-14
Status: phase 1 (local frames, `elliptical_disc` and `ellipse` patches and carveouts, prism carveouts, rotation parsing) implemented on this branch; 3D (phase 2) still a proposal

## Summary

Add elliptical shapes as new patch and carveout types, named and parameterized exactly as in PhysiCellECMCreator.jl (`elliptical_disc`, `ellipse` with `thickness`, `a`, `b`, `rotation`). Do not overload `disc` and `annulus` with axis parameters; keep them as circle shorthands that resolve to the same internal geometry.

Underneath, restructure placement around one idea: **every shape lives in a local frame** given by a center, a `normal`, and a `rotation` about that normal. Shape parameters are in local coordinates. Sampling, lattice generation, and inside tests are written once per shape in local coordinates, and the frame maps them into the domain. That one rule is what makes the 3D step later (spheres, ellipsoids, boxes, volumetric `everywhere`) an addition rather than a rewrite.

Sequencing: land the hexagonal grid PR first, then ellipses (phase 1), then 3D (phase 2). Details and the decisions you need to make are below.

## Where things stand

### This package (main, v0.0.9)

| Patch type | Parameters | Orientation | Carveout version |
|---|---|---|---|
| `disc` | `x0 y0 z0 radius` | `normal` (any plane) | yes, but tests a **ball** (3D distance), not a disc |
| `annulus` | `x0 y0 z0 inner_radius outer_radius` | `normal` | yes, tests a spherical shell |
| `rectangle` | corner `x0 y0 z0` + `width height` or `x1 y1` | xy-plane only | yes, tests an infinite prism in z |
| `everywhere` | `number` | xy-plane at z = 0, even in a 3D domain | no |

Observations that shape the proposal:

- Disc and annulus share one sampler (`placeAnnulus!`) with a radius function. Rectangle and everywhere each have their own. There is no shape abstraction; carveouts are a parallel set of structs with their own inside tests.
- The carveout semantics are inconsistent: circular carveouts are balls, rectangular carveouts are prisms. For the common case (everything in the z = z0 plane) they agree, so nobody has noticed.
- On main, a disc with no `normal` and `z0 != 0` lands at `2 * z0` (z0 is added twice). The hexagonal grid PR fixes this.
- The in-plane basis for a tilted disc is discontinuous at the default normal: for `normal = (0, 0, 1)` the basis is the x and y axes, but for `normal = (ε, 0, 1)` the first basis vector is `(0, -1, 0)`. Circles do not care. Ellipses with a `rotation` will, so the frame construction has to change before ellipses can be added.
- `parseCarveoutPatchCollection` hardcodes the carveout type list instead of calling `supportedCarveoutTypes()`. Every new type has to be added in two places.

### The open hexagonal grid PR (drbergman-lab/PhysiCellCellCreator.jl#31)

It touches the same code this brief would touch. It adds a `placement` attribute (`random` or `hexagonal`), `cell_radius` and `spacing_factor` keywords on `generateICCell`, and refactors placement into `annulusCoordinates`, `placeBox!`, `filterCells!`, and `hexagonalGridPoints`. Its lattice is generated in local `(c1, c2)` coordinates and then mapped by the normal, which is already the local-frame idea in embryo.

That PR merged on 2026-10-14 and phase 1 was built on top of it.

### PhysiCellECMCreator.jl conventions to match

- Types: `everywhere`, `ellipse` (a ring of `thickness` extending **outward** from the ellipse `(a, b)` to `(a + thickness, b + thickness)`), `elliptical_disc` (filled), `ellipse_with_shell` (interior + shell + optional exterior in one patch).
- Elements: `x0 y0 a b rotation thickness`. `a` is along local x, `b` along local y.
- `rotation` is counter-clockwise; degrees by default; radians if the text contains `pi` or `π`; a `units` attribute starting with `rad` forces radians. Only 2D: no normal.
- Internally an `EllipticalECM` struct holds the center, axes, rotation, thickness, and the two rotation matrices, with an `_ellipseFunctional(x, y, a, b) = (x/a)^2 + (y/b)^2` predicate and `positionRelativeToEllipse` returning `:inside`, `:thickness`, or `:outside`. This is the right shape of API to reuse.
- One thing not to copy: `parseRotation` evaluates `pi/6` with `eval(Meta.parse(...))`. Replace that with a small parser for `[k*]pi[/m]` forms (with optional sign). Same accepted inputs, no code execution from an XML file.

### How PCMM calls these

PCMM reads `x/y/z_min/max` from the config and passes all six to `generateICCell`, so the 3D domain is already available. For ECM it passes only x and y plus `dx`, `dy`, with a comment that 3D is not yet supported. PCMM does not currently pass PhysiCell's `use_2D` flag to either package.

## Proposal

### 1. One internal model: shape in a local frame

```julia
struct Frame
    center::Vector{Float64}   # (x0, y0, z0)
    R::Matrix{Float64}        # 3x3, columns are local axes u1, u2, n
end

abstract type Shape end
abstract type PlanarShape <: Shape end      # 2D: local coords (c1, c2), embedded at local c3 = 0
abstract type SolidShape <: Shape end       # 3D: local coords (c1, c2, c3)

struct EllipticalDisc <: PlanarShape; a; b; end
struct EllipticalRing <: PlanarShape; a; b; thickness; end   # ECM's `ellipse`
struct Rectangle      <: PlanarShape; width; height; end     # corner-at-origin, see decision D4
```

Each shape implements three methods in local coordinates:

| Method | Used by |
|---|---|
| `sampleUniform(shape, n)` | random placement |
| `inside(shape, c)` | carveouts, lattice filtering, and shells done by rejection |
| `bounds(shape)` | hexagonal (later HCP) lattice generation |

The frame does the rest: `toDomain(frame, local_points)` and `toLocal(frame, point)`. Carveouts are `(Shape, Frame)` pairs whose test is `inside(shape, toLocal(frame, p))`. A `PlanarShape` used as a carveout ignores the local third coordinate, i.e. it is an infinite prism along its normal (decision D3).

Existing types map onto this with no new behavior: `disc` is `EllipticalDisc(r, r)`, `annulus` is `EllipticalRing(r_in, r_in, r_out - r_in)`, `rectangle` is `Rectangle(w, h)` with the frame centered at the corner.

### 2. Frame construction: `normal` and `rotation`

Given unit normal `n` and rotation angle `θ`:

1. `u1` is the unit projection of the x-axis onto the plane: `normalize(e_x - (e_x · n) n)`. If `n` is parallel to the x-axis, use the y-axis instead.
2. `u2 = n × u1`.
3. Rotate `u1`, `u2` counter-clockwise about `n` by `θ`.

For the default normal this gives `u1 = e_x`, `u2 = e_y`, so `rotation` means exactly what it means in ECM Creator, and a slightly tilted normal gives a slightly tilted frame. The guide should state the rule in one sentence: "`a` runs along the x-axis projected onto the patch plane, `b` is perpendicular to it in the plane, and `rotation` turns both counter-clockwise about the normal."

This replaces the basis in `annulusCoordinates` (hexagonal PR) for every shape, including the hexagonal lattice orientation.

### 3. New XML types (phase 1)

Patch types, matching ECM Creator names and elements, plus the elements this package already uses on every patch:

```xml
<patch_collection type="elliptical_disc">
  <patch ID="1">
    <x0>0.0</x0> <y0>0.0</y0> <z0>0.0</z0>
    <a>120.0</a> <b>60.0</b>
    <rotation>pi/8</rotation>          <!-- same parsing rules as ECM Creator -->
    <normal>0,0,1</normal>             <!-- optional, as for disc -->
    <number>200</number>               <!-- or placement="hexagonal" + spacing -->
  </patch>
</patch_collection>

<patch_collection type="ellipse">
  <patch ID="1">
    <x0>0.0</x0> <y0>0.0</y0> <z0>0.0</z0>
    <a>120.0</a> <b>60.0</b>
    <thickness>30.0</thickness>        <!-- ring from (a,b) out to (a+30, b+30) -->
    <rotation units="deg">30</rotation>
    <number>80</number>
  </patch>
</patch_collection>
```

Carveout types `elliptical_disc` and `ellipse` with the same geometric elements and no `number`. This is the "elliptical carve out" you asked for, e.g. a tumor-shaped hole in a fibroblast field.

`ellipse_with_shell` is deliberately **not** proposed. In ECM Creator it exists because one patch has to assign three regions with different properties. Here a patch places one cell type, so the same thing is an `elliptical_disc` patch for the interior plus an `ellipse` patch (same `x0 y0 a b rotation`) for the shell, possibly under different `cell_patches`. Say so in the guide so the two packages still feel parallel.

`disc` and `annulus` stay as they are. Do not add `a`/`b` to them. Two reasons: a "disc" with unequal axes is misleading, and ECM Creator already picked the names, so `elliptical_disc` is what users of both packages will look for.

### 4. Sampling

- **Elliptical disc**: sample the unit disc (`r = sqrt(u)`, uniform `θ`), scale by `(a, b)`, map through the frame. Exact and uniform, since scaling is affine.
- **Ring with `thickness`**: the ECM convention `(a + t, b + t)` is an offset, not a scaling, so there is no affine trick. Sample the outer elliptical disc and reject points inside the inner ellipse; the existing `createCellsDataFrame` retry loop with `max_fails` already handles this (it re-samples only the shortfall). Expected acceptance is `1 - ab / ((a+t)(b+t))`, which is fine even for thin rings because each round is vectorized.
- **Special case worth keeping**: when the inner ellipse is a scaling of the outer (`a == b`, i.e. today's `annulus`), use the exact radial sampler `r = sqrt(r_in^2 + (r_out^2 - r_in^2) u)` so `annulus` behavior is unchanged bit for bit (modulo RNG stream).
- **Hexagonal placement**: generate the lattice over `bounds(shape)` in local coordinates, keep points with `inside(shape, c)`, map through the frame. The hexagonal PR already does this for circles with an explicit radius test; `inside` generalizes it.

### 5. Code organization

The single source file is about 400 lines today and roughly 800 after the hexagonal PR. Adding shapes is the moment to split it:

```
src/PhysiCellCellCreator.jl   # module, exports, generateICCell
src/frame.jl                  # Frame, normal/rotation parsing, to/from local
src/shapes.jl                 # Shape types, sampleUniform, inside, bounds
src/parse.jl                  # XML element parsing (center, rotation, thickness, bools, number)
src/patches.jl                # patch type registry, generatePatch!, placement (random/hexagonal)
src/carveouts.jl              # carveout registry (reuses shapes + frames)
src/template.jl               # createICCellXMLTemplate
```

Both patch and carveout registries should be driven from one table of `name => (ShapeConstructor, required elements)` so `supportedPatchTypes`, `supportedCarveoutTypes`, and the two `supported*TextElements` functions are derived rather than maintained by hand, and the hardcoded list in `parseCarveoutPatchCollection` goes away.

### 6. Tests to add

- Uniformity: for an `elliptical_disc` with 10^5 cells, the fraction inside the ellipse scaled by `s` is `s^2 ± tolerance`; for an `ellipse` ring, the fraction inside the inner ellipse is 0.
- Orientation: with `a >> b` and `rotation = θ`, the principal axis of the sample covariance is within a degree of `θ`; repeat with a tilted `normal` and check the samples lie in the plane (`|(p - center) · n| < 1e-10`).
- Frame continuity: `normal = (1e-6, 0, 1)` and `normal = (0, 0, 1)` give the same principal axis.
- Carveouts: no cell inside an `elliptical_disc` or `ellipse` carveout; a rotated ellipse carveout removes the right fraction from a rectangle patch.
- Hexagonal: lattice count inside an ellipse is within a few percent of `π a b / (spacing^2 √3 / 2)`.
- Rotation parsing: `30`, `pi/6`, `π/6`, `units="rad"`, `units="deg"`, negative values, and the rejected combination `pi` with `units="deg"`, without `eval`.
- Regression: template XML still generates; `disc`/`annulus`/`rectangle` samples unchanged in distribution.

## 3D roadmap (phase 2)

Everything above is chosen so that 3D is additive. What changes:

**Dimensionality signal.** `domain_dict` gains an optional `"use_2D"` key (same name as the PhysiCell config field, so PCMM can pass it straight through), defaulting to `true` for backward compatibility. It affects only `everywhere`: in 3D it fills the full box `[x_min, x_max] × [y_min, y_max] × [z_min, z_max]` instead of the z = 0 plane. Every other shape already carries its own dimensionality in its type. PCMM needs a one-line change to read `use_2D` from the config.

**Solid shapes** as `SolidShape` subtypes, each with the same three methods:

| Type | Local parameters | Sampler |
|---|---|---|
| `sphere` | `radius` | uniform direction × `r = radius · u^(1/3)` |
| `spherical_shell` | `inner_radius outer_radius` | exact: `r = (r_in^3 + (r_out^3 - r_in^3) u)^(1/3)` |
| `ellipsoid` | `a b c` | unit ball scaled by `(a, b, c)` |
| `ellipsoidal_shell` | `a b c thickness` | rejection from the outer ellipsoid (offset, as for `ellipse`) |
| `box` | corner + `width height depth` (or `x1 y1 z1`) | uniform in each axis |

All take `normal` and `rotation` through the same `Frame`: `normal` is the local c-axis, `rotation` spins `a`/`b` about it. That is the whole orientation story, consistent with the planar shapes.

Each one is also a carveout, and a solid carveout tests the full local coordinate, which is how `disc`-carveout-as-ball gets an explicit name (`sphere`) instead of being an accident.

**Naming note.** ECM Creator's convention (curve name = ring, `*_disc` = filled) does not extend to 3D, where `sphere` and `ellipsoid` are universally read as solids in modeling contexts. The table above uses the solid reading plus an explicit `_shell` suffix. Flag this to ECM Creator when it goes 3D so the two packages pick the same 3D vocabulary.

**Hexagonal → close packing.** The `placement` attribute from the hexagonal PR is already the extension point: a later `placement="hcp"` (or `"fcc"`) generates a 3D lattice over `bounds(shape)` and filters with `inside`. Not part of this proposal, just confirming nothing here blocks it.

**ECM Creator 3D** is a separate, larger change (voxel grid in z, orientation vectors in 3D) and is not in scope. Sharing the shape vocabulary is the only coupling.

## Decisions for you

- **D1. Ring parameterization.** Recommended: `thickness` extending outward, exactly as ECM Creator. Alternative: `inner_a inner_b outer_a outer_b`, more general but diverges from ECM and invites inconsistent aspect ratios.
- **D2. Keep `disc`/`annulus` as separate types** (recommended) vs. folding them into `elliptical_disc`/`ellipse` with `radius` as shorthand. Separate types cost nothing once both resolve to the same shape.
- **D3. Carveout semantics for planar shapes.** Recommended: infinite prism along the carveout's normal (matches `rectangle` today, and solids are what `sphere` carveouts are for). This changes `disc`/`annulus` carveouts from balls to cylinders for cells not in the carveout's plane. Nobody is likely relying on that, but it goes in the changelog. Alternative: keep the ball behavior and document it as a quirk.
- **D4. Rectangles and frames.** `rectangle` is corner-based while every other shape is center-based. Recommended for phase 1: leave `rectangle` alone (no `normal`/`rotation`). If rotation is added later, rotate about the rectangle's center, not the corner, and say so. Alternative: switch rectangles to center-based now with `x1`/`y1` kept as the corner form, which is a breaking change in the XML.
- **D5. Shared parsing with ECM Creator.** Recommended: copy the rotation/thickness parsing (minus `eval`) into this package for now. Extracting a tiny shared package (say `PhysiCellICBase.jl`) is only worth it once a third consumer exists or the two copies drift.
- **D6. Sequencing.** Recommended: merge the hexagonal grid PR, then phase 1 in one PR (frame refactor, shapes, elliptical patches and carveouts, template, guide, tests), then 3D in a second PR that also carries the PCMM `use_2D` change.

## Out of scope, noted for later

- Cells on a curve (`thickness = 0` on an `ellipse` meaning "on the perimeter"), which needs arc-length sampling.
- Non-uniform densities inside a patch (radial gradients).
- Polygon / arbitrary mask patches from an image.
