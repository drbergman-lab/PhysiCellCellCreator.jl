```@meta
CurrentModule = PhysiCellCellCreator
```

# Adding a new shape

This page is for developers adding a new patch (and carveout) type.
Placement is built around two ideas, and a new shape only has to fit them:

- A **[`Shape`](@ref)** is a region defined in its own local coordinates `(c1, c2)` (planar shapes live in the plane `c3 = 0`).
- A **[`Frame`](@ref)** positions the shape in the domain: a center and an orthonormal basis `(u1, u2, n)`.

Random placement, hexagonal placement, carveouts, and `restrict_to_domain` are written once in terms of three shape methods and the frame, so a new shape gets all of them for free.
The `line` patch (`Segment` in `src/shapes.jl`, `LinePatch` in `src/patches.jl`) is a small, complete example of every step below.

## 1. Define the shape (`src/shapes.jl`)

Add a struct that subtypes `PlanarShape` and validates its parameters in an inner constructor:
```julia
struct Segment <: PlanarShape
    length::Float64
    function Segment(length::Real)
        @assert length > 0 "The length of a line must be positive. Got $(length)."
        return new(length)
    end
end
```

Then implement the three methods every shape needs, all in local coordinates:

| Method | Returns | Used by |
|---|---|---|
| `sampleUniform(shape, n)` | `(c1, c2)` vectors of `n` points uniform over the shape | random placement |
| `inside(shape, c1, c2; tol=0.0)` | `BitVector` of which points are in the shape (boundary included) | carveouts, lattice filtering, rejection sampling |
| `bounds(shape)` | `(half_width, half_height)` of the local bounding box, centered on the origin | hexagonal lattice generation |

Rules for these:
- The shape is centered on its local origin: hexagonal placement puts a lattice point there, and `bounds` is symmetric about it.
- `inside` is the **only** membership test in the package. Do not write a second one elsewhere.
- `tol` enlarges the shape by that distance (for a ring, it also shrinks the hole). It is passed only when filtering lattice points, to keep points that land on the boundary up to floating-point rounding. Random sampling and carveouts call `inside` with `tol = 0`.
- If a shape has to sample by rejection (as `EllipticalRing` does), reject with `inside` and cap each proposal batch with `max_proposal_batch`.

## 2. Register the patch type (`src/patches.jl`)

Three additions:

1. A marker type and a registry entry. The registry order is the order `supportedPatchTypes()` reports:
   ```julia
   struct LinePatch end
   const patch_type_registry = [
       ...,
       "line" => LinePatch,
   ]
   ```
2. The geometric XML elements of the type. These are also the elements of the carveout of the same name, and `supportedPatchTextElements` appends the placement elements (`number`, `spacing`, `max_fails`, `restrict_to_domain`) automatically:
   ```julia
   geometryElements(::Type{LinePatch}) = ["x1", "y1", "z1", "x2", "y2", "z2"]
   ```
   PhysiCellModelManager.jl validates variation targets against these lists, so list exactly the elements the parser reads.
3. `parseGeometry`, which turns the XML into a `(shape, frame)` pair:
   ```julia
   function parseGeometry(::Type{LinePatch}, patch::XMLElement, descriptor::String)
       p1 = [parseFloatElement(patch, "x1", descriptor), parseFloatElement(patch, "y1", descriptor), parseOptionalFloatElement(patch, "z1", 0.0)]
       p2 = [parseFloatElement(patch, "x2", descriptor), parseFloatElement(patch, "y2", descriptor), parseOptionalFloatElement(patch, "z2", 0.0)]
       @assert p1 != p2 "The endpoints of line $(descriptor) coincide ($(p1)). A line must have positive length."
       return Segment(norm(p2 .- p1)), frameAlong((p1 .+ p2) ./ 2, p2 .- p1)
   end
   ```
   Use the helpers in `src/parse.jl` (`parseFloatElement`, `parsePositiveFloatElement`, `parseOptionalFloatElement`, `parseCenter`, ...) so that missing or invalid elements produce an `AssertionError` naming the patch (`descriptor` is "cell type with patch ID n", or "carveout ID m of ..." for carveouts).

The frame comes from one of:
- `parseCenteredFrame(patch, descriptor)`: center `x0 y0 z0` plus optional `normal` and `rotation`. Use this for any shape with a center and an orientation; it is what the disc, annulus, and elliptical patches use, and it gives the shape the [Orientation](@ref) rules documented for users.
- `frameAlong(center, direction)`: local x along a direction, for shapes that only use their local x axis (the line).
- `Frame(center, normal, rotation)` directly, for a fixed orientation (the rectangle and `everywhere` use the default normal and zero rotation).

The generic `generatePatch!` then handles placement. Only `everywhere` has its own method, because its geometry comes from the domain rather than the XML.

## 3. Carveouts (`src/carveouts.jl`)

Nothing to add: every registered patch type is a carveout type with the same geometric elements, and a planar carveout removes cells whose projection onto its plane is `inside` the shape (an infinite prism along its normal).
If the new type must not be a carveout, add it to the exclusion list in `supportedCarveoutTypes()` with a reason (`everywhere` would remove every cell; `line` has no area and would remove none).

## 4. Template, guide, and tests

- Add an example of the new type to `createICCellXMLTemplate` in `src/template.jl` (the `addPatchCollection!`, `addPatch!`, and `addCarveout!` helpers keep this short). The template is what `PhysiCellModelManager.createICCellXMLTemplate` hands to users.
- Document the type in the [Guide](@ref): its elements, an XML example, how hexagonal placement behaves on it, and whether it can be a carveout. Update the Example XML at the end to match the template.
- Add a test file `test/<Name>Tests.jl` and include it from `test/runtests.jl`. The tests for existing shapes are the pattern; cover at least:
  - random placement: every cell is `inside` the shape, and a sub-region holds the expected fraction of cells (uniformity);
  - orientation, if the shape has one: `normal` and `rotation` place it where the guide says;
  - hexagonal placement: a cell at the center, nearest-neighbor distance equal to `spacing`, cell count within a few percent of area divided by lattice cell area, boundary points kept;
  - the shape as a carveout: no surviving cell inside it, and the right number removed from a lattice;
  - `restrict_to_domain` and carveouts applied to the new patch;
  - invalid XML (missing or non-positive elements) throws an `AssertionError`;
  - the registries: `supportedPatchTypes()`, `supportedPatchTextElements("<name>")`, and the carveout lists.

## Beyond planar shapes

Solid (3D) shapes are not implemented yet. The intended route is a `SolidShape <: Shape` whose methods take `(c1, c2, c3)`: `sampleUniform` returns three coordinates, `bounds` returns three half-extents, `inside` tests the full local point (so a solid carveout is a true volume, not a prism), and `placeShape!` gets a method for `SolidShape` that maps all three coordinates through the frame.
`Frame` already carries a full 3D basis, so no change is needed there: `normal` is the local c3 axis and `rotation` spins the other two about it.
