# Guide

You can create a simple XML by running
```julia
using PhysiCellCellCreator
path_to_folder = "template_folder"
createICCellXMLTemplate(path_to_folder) # will create cells.xml inside the folder "template_folder"
```

Then, you can edit the XML file as needed.
Note, the root of the XML file is `<ic_cells>`.

## Cell patches
The children of the root are `<cell_patches>` elements, each corresponding to a different cell type in the PhysiCell simulation.
The cell type is set by the attribute `name`, e.g., `<cell_patches name="fibroblast">`.

## Patch collection
Each `<cell_patches>` element contains a collection of `<patch_collection>` elements, each containing a particular type of patch.
The currently supported patch types are `everywhere`, `disc`, `annulus`, `elliptical_disc`, `ellipse`, `rectangle`, and `line`.
A typical example is `<patch_collection type="everywhere">`.

The elliptical patches use the same names, elements, and conventions as the patches of [PhysiCellECMCreator.jl](https://github.com/drbergman-lab/PhysiCellECMCreator.jl), so an ellipse described in one package describes the same region in the other.

## Patch
Each `<patch_collection>` element contains a collection of `<patch>` elements.
Each patch element has an `ID` attribute, which can be used to identify the patch.
Depending on the type of patch, the patch may have different parameters.

### `everywhere`
The `everywhere` patch just has a `number` element, specifying the number of cells to place everywhere.
```xml
<patch_collection type="everywhere">
  <patch ID="1">
      <number>100</number>
  </patch>
</patch_collection>
```
Cells are placed in the `xy`-plane at `z = 0`, even if the domain is 3D.

### `disc`
The `disc` patch has the following parameters:
- `x0`, `y0`, `z0`: the center of the disc
- `radius`: the radius of the disc
- `number`: the number of cells to place in the disc
```xml
<patch_collection type="disc">
  <patch ID="1">
      <x0>100.0</x0>
      <y0>100.0</y0>
      <z0>0.0</z0>
      <radius>50.0</radius>
      <number>50</number>
  </patch>
</patch_collection>
```

A `disc` patch can also have `normal` and `rotation` elements; see [Orientation](@ref).
Note: the `disc` really is a 2D disc, not a solid sphere even if the PhysiCell simulation is 3D.

### `annulus`
The `annulus` patch has identical parameters to the `disc` patch, but with `radius` replaced by `inner_radius` and `outer_radius`.
```xml
<patch_collection type="annulus">
  <patch ID="1">
      <x0>100.0</x0>
      <y0>100.0</y0>
      <z0>0.0</z0>
      <inner_radius>50.0</inner_radius>
      <outer_radius>100.0</outer_radius>
      <number>50</number>
  </patch>
</patch_collection>
```

### `elliptical_disc`
The `elliptical_disc` patch is a filled ellipse with the following parameters:
- `x0`, `y0`, `z0`: the center of the ellipse
- `a`, `b`: the semi-axes of the ellipse along the patch's local `x` and `y` directions, respectively (see [Orientation](@ref))
- `rotation`: optional, the counter-clockwise rotation of the ellipse about its normal (see [Orientation](@ref))
- `number`: the number of cells to place in the ellipse

With the default orientation, the patch is the region

```math
\begin{aligned}
\frac{(x-x_0)^2}{a^2} + \frac{(y-y_0)^2}{b^2} \le 1
\end{aligned}
```

rotated counter-clockwise by `rotation` about its center.
```xml
<patch_collection type="elliptical_disc">
  <patch ID="1">
      <x0>100.0</x0>
      <y0>100.0</y0>
      <z0>0.0</z0>
      <a>80.0</a>
      <b>40.0</b>
      <rotation>pi/6</rotation>
      <number>50</number>
  </patch>
</patch_collection>
```
A `disc` is an `elliptical_disc` with `a` and `b` both equal to `radius`.

### `ellipse`
The `ellipse` patch is a ring of cells around an ellipse: it has the same parameters as `elliptical_disc` plus a `thickness`.
Cells are placed between the ellipse with semi-axes `a` and `b` and the ellipse with semi-axes `a + thickness` and `b + thickness`, i.e., the ring extends outward from the ellipse.
With the default orientation, the patch is the region

```math
\begin{aligned}
1 \le \frac{(x-x_0)^2}{a^2} + \frac{(y-y_0)^2}{b^2}
\quad\text{and}\quad
\frac{(x-x_0)^2}{(a+\text{thickness})^2} + \frac{(y-y_0)^2}{(b+\text{thickness})^2} \le 1
\end{aligned}
```

```xml
<patch_collection type="ellipse">
  <patch ID="1">
      <x0>100.0</x0>
      <y0>100.0</y0>
      <z0>0.0</z0>
      <a>80.0</a>
      <b>40.0</b>
      <thickness>20.0</thickness>
      <rotation>pi/6</rotation>
      <number>50</number>
  </patch>
</patch_collection>
```
An `annulus` is an `ellipse` with `a` and `b` both equal to `inner_radius` and `thickness` equal to `outer_radius - inner_radius`.

To fill an ellipse with one cell type and surround it with a shell of another (the `ellipse_with_shell` patch of PhysiCellECMCreator.jl), use an `elliptical_disc` patch under one `cell_patches` element and an `ellipse` patch with the same `x0`, `y0`, `z0`, `a`, `b`, and `rotation` under the other.

### `rectangle`
The `rectangle` patch is assumed to be parallel to the xy-plane and has the following parameters:
- `x0`, `y0`: the bottom-left corner of the rectangle (minimum x and y)
- `z0`: the z-coordinate of the rectangle
- `width`, `height`: the width and height of the rectangle in the x and y directions, respectively
- `number`: the number of cells to place in the rectangle
```xml
<patch_collection type="rectangle">
  <patch ID="1">
      <x0>100.0</x0>
      <y0>100.0</y0>
      <z0>0.0</z0>
      <width>50.0</width>
      <height>100.0</height>
      <number>50</number>
  </patch>
</patch_collection>
```

Note: this patch will place cells with x-coordinate in [100.0, 150.0] and y-coordinate in [100.0, 200.0].
Instead of `width` and `height`, the opposite corner can be given with `x1` and `y1`.

### `line`
The `line` patch places cells on the straight line segment between two points:
- `x1`, `y1`, `z1`: the first endpoint
- `x2`, `y2`, `z2`: the second endpoint
- `number`: the number of cells to place on the line

`z1` and `z2` are optional and default to `0`, so a 2D line needs only `x1`, `y1`, `x2`, and `y2`.
```xml
<patch_collection type="line">
  <patch ID="1">
      <x1>-100.0</x1>
      <y1>0.0</y1>
      <x2>100.0</x2>
      <y2>50.0</y2>
      <number>20</number>
  </patch>
</patch_collection>
```
With random placement, cells are placed uniformly along the segment.
With `placement="hexagonal"` (see [Hexagonal grid](@ref)), cells are evenly spaced along the segment at the given `spacing`, with a cell at the midpoint; the endpoints get cells only if the half-length is a multiple of the spacing.
A `line` has no `normal` or `rotation`: its direction is set by its endpoints.
Carveouts apply as usual, so a line passing through a carveout gets a gap.

## Orientation
The `disc`, `annulus`, `elliptical_disc`, and `ellipse` patches (and the carveouts of the same names) can be oriented in 3D with two optional elements.

The `normal` element sets the normal vector of the plane the patch lies in, as three comma-separated numbers, e.g., `<normal>1.0,0.0,0.0</normal>`.
If omitted, the normal vector is assumed to be `(0, 0, 1)`, i.e., the patch lies in the `xy`-plane.
The vector does not need to be normalized.

Within its plane, the patch's local `x` direction (along which `a` is measured) is the projection of the domain's `x`-axis onto the plane, or of the `y`-axis if the normal is parallel to the `x`-axis.
The local `y` direction (along which `b` is measured) is perpendicular to it in the plane, such that the local `x`, local `y`, and normal directions form a right-handed system.
For the default normal, these are just the domain's `x` and `y` axes.

The `rotation` element then turns the patch counter-clockwise about its normal, i.e., counter-clockwise as seen by a viewer at the tip of the normal looking back at the patch.
For the default normal this is counter-clockwise as seen from above (from `+z`); for the normal `(0, 0, -1)` the same `rotation` value turns the patch clockwise as seen from above.
Its value follows the rules of PhysiCellECMCreator.jl:
- by default, the value is in degrees, e.g., `<rotation>30</rotation>`;
- if the value contains `pi` or `π`, it is in radians, e.g., `<rotation>pi/6</rotation>`, `<rotation>-2pi/3</rotation>`, or `<rotation>0.5*π</rotation>`;
- the units can be set explicitly with a `units` attribute starting with `rad` or `deg`, e.g., `<rotation units="rad">0.5236</rotation>`.

Using `pi` or `π` with `units="deg"` is an error.
For a `disc` or `annulus`, the rotation only matters for a [Hexagonal grid](@ref), whose lattice turns with the patch.

## Refined control of cell placement

### Restrict to domain
Each of the above patches (except `everywhere`) can have a `restrict_to_domain` element, which will force cells to be placed within the domain as specified by the `domain_dict` parameter in the `generateICCell` function.
If omitted, it is assumed to be "true".
Possible values are "true", "false", "0", and "1".
Add the following to the patch element to set: `<restrict_to_domain>false</restrict_to_domain>`.

### Hexagonal grid
Each `<patch>` element can carry a `placement` attribute saying how cells are placed in it.
The default, `placement="random"`, places `number` cells uniformly at random as described above.
Setting `placement="hexagonal"` instead fills the patch with a hexagonal grid (lattice) of cells; the attribute tells you what to look for inside the patch: a `spacing` element rather than a `number` element.
The grid is centered on the patch (the center of a `disc`, `annulus`, `elliptical_disc`, or `ellipse`, the center of a `rectangle`, the midpoint of a `line`, or the `x`/`y` center of the domain at `z = 0` for `everywhere`, as with random `everywhere` placement) with a lattice point at that center.
For a `line`, the grid reduces to evenly spaced cells along the line.
Only the lattice points inside the patch are kept, and `restrict_to_domain` and carveouts are applied as usual.
The grid lies in the plane of the patch and its rows run along the patch's local `x` direction, so `normal` and `rotation` orient the lattice along with the patch (see [Orientation](@ref)).

The distance between neighboring cell centers is set by the optional `spacing` element:
```xml
<patch_collection type="disc">
  <patch ID="1" placement="hexagonal">
      <x0>100.0</x0>
      <y0>100.0</y0>
      <z0>0.0</z0>
      <radius>50.0</radius>
      <spacing>12.0</spacing>
  </patch>
</patch_collection>
```

If `spacing` is omitted, the spacing is computed from the cell radius passed to `generateICCell` via the `cell_radius` keyword argument as `spacing = spacing_factor * 2 * cell_radius`.
This lets a calling package that knows the cell radius (e.g., from the PhysiCell configuration file) fill patches without hard-coding the spacing in the XML.
`cell_radius` can be a single number applied to all cell types or a `Dict` mapping cell type names to radii.
The `spacing_factor` keyword argument defaults to `0.95`, matching the slightly compressed packing used in PhysiCell's sample projects; use `1.0` to place cells tangent to their neighbors.
```julia
generateICCell(path_to_ic_cell_xml, domain_dict; cell_radius=8.4)
generateICCell(path_to_ic_cell_xml, domain_dict; cell_radius=Dict("default" => 8.4, "fibroblast" => 6.2))
generateICCell(path_to_ic_cell_xml, domain_dict; cell_radius=8.4, spacing_factor=1.0)
```
An error is thrown if a hexagonal grid patch has no `spacing` and no radius is available for its cell type.

Notes:
- A `placement="hexagonal"` patch cannot have a `number` element as both determine how many cells are placed. The number of cells in a hexagonal grid is determined by the patch size and the spacing.
- A random placement patch cannot have a `spacing` element.
- `max_fails` is ignored for hexagonal grids since placement is deterministic. A warning is issued if no lattice points survive the patch, domain, and carveout filters.
- Lattice points exactly on the patch boundary are kept. To avoid floating-point rounding dropping them, points within `boundary_tolerance * spacing` outside the boundary are also kept; `boundary_tolerance` is a keyword argument of `generateICCell` defaulting to `1e-10`.

### Carveouts
Every patch type can also include a `<carveout_patches>` element, which specifies regions where cells should not be placed.
The currently supported carveouts are `disc`, `annulus`, `elliptical_disc`, `ellipse`, and `rectangle`; `everywhere` and `line` cannot be carveouts (the first would remove every cell, the second none).
These follow the same structure as the patches above, but do not contain a `number` (or `spacing`) element.

A carveout removes every cell whose projection onto the carveout's plane lies inside the carveout's shape, i.e., a carveout is an infinite prism along its normal.
`disc`, `annulus`, `elliptical_disc`, and `ellipse` carveouts may have `normal` and `rotation` elements, as described in [Orientation](@ref); by default they are prisms along the `z`-axis.
`rectangle` carveouts are always prisms along the `z`-axis.
For cells in the plane of the carveout (the usual case), this is just "inside the shape".

To add a carveout to a patch, add a `<carveout_patches>` element to the patch, then follow the same XML structure as the patches above.
```xml
<patch_collection type="disc">
  <patch ID="1">
      <x0>100.0</x0>
      <y0>100.0</y0>
      <z0>0.0</z0>
      <radius>50.0</radius>
      <number>50</number>
      <carveout_patches>
          <patch_collection type="elliptical_disc">
              <patch ID="1">
                  <x0>130.0</x0>
                  <y0>130.0</y0>
                  <z0>0.0</z0>
                  <a>15.0</a>
                  <b>8.0</b>
                  <rotation>45</rotation>
              </patch>
          </patch_collection>
      </carveout_patches>
  </patch>
</patch_collection>
```

### Max fails
When attempting to place cells, points taken uniformly randomly from the entire patch are selected and then pared down based on the domain and carveouts.
This can lead to a situation where many samplings must be taken to finally achieve `number` of cells in each patch.
In the worst case scenario, a patch has no area and cannot support any cells.
To prevent the program from running indefinitely, a `max_fails` parameter controls how many samples are taken before giving up.
If omitted, it is assumed to be 100.
Add the following to the patch element (not the carveout, but the patch itself) to set: `<max_fails>100</max_fails>`.

## Example XML
A full example XML, as created by `createICCellXMLTemplate`, is shown below.
```xml
<?xml version="1.0" encoding="utf-8"?>
<ic_cells>
  <cell_patches name="default">
    <patch_collection type="disc">
      <patch ID="1">
        <x0>0.0</x0>
        <y0>0.0</y0>
        <z0>0.0</z0>
        <radius>40.0</radius>
        <number>50</number>
        <normal>0,0,1</normal>
        <max_fails>100</max_fails>
        <carveout_patches>
          <patch_collection type="rectangle">
            <patch ID="1">
              <x0>0.0</x0>
              <y0>0.0</y0>
              <z0>0.0</z0>
              <width>10.0</width>
              <height>10.0</height>
            </patch>
          </patch_collection>
        </carveout_patches>
      </patch>
      <patch ID="2" placement="hexagonal">
        <x0>200.0</x0>
        <y0>200.0</y0>
        <z0>0.0</z0>
        <radius>40.0</radius>
        <spacing>12.0</spacing>
      </patch>
    </patch_collection>
    <patch_collection type="annulus">
      <patch ID="1">
        <x0>50.0</x0>
        <y0>50.0</y0>
        <z0>0.0</z0>
        <inner_radius>10.0</inner_radius>
        <outer_radius>200.0</outer_radius>
        <number>50</number>
        <restrict_to_domain>true</restrict_to_domain>
      </patch>
    </patch_collection>
    <patch_collection type="rectangle">
      <patch ID="1">
        <x0>-50.0</x0>
        <y0>-50.0</y0>
        <z0>0.0</z0>
        <x1>50.0</x1>
        <y1>50.0</y1>
        <number>10</number>
        <carveout_patches>
          <patch_collection type="disc">
            <patch ID="1">
              <x0>0.0</x0>
              <y0>0.0</y0>
              <z0>0.0</z0>
              <radius>10.0</radius>
            </patch>
          </patch_collection>
        </carveout_patches>
      </patch>
    </patch_collection>
    <patch_collection type="elliptical_disc">
      <patch ID="1">
        <x0>-200.0</x0>
        <y0>200.0</y0>
        <z0>0.0</z0>
        <a>80.0</a>
        <b>40.0</b>
        <rotation>pi/6</rotation>
        <number>60</number>
      </patch>
    </patch_collection>
    <patch_collection type="ellipse">
      <patch ID="1" placement="hexagonal">
        <x0>-200.0</x0>
        <y0>200.0</y0>
        <z0>0.0</z0>
        <a>80.0</a>
        <b>40.0</b>
        <thickness>20.0</thickness>
        <rotation>pi/6</rotation>
        <spacing>10.0</spacing>
      </patch>
    </patch_collection>
    <patch_collection type="line">
      <patch ID="1" placement="hexagonal">
        <x1>-400.0</x1>
        <y1>-450.0</y1>
        <x2>400.0</x2>
        <y2>-350.0</y2>
        <spacing>20.0</spacing>
      </patch>
    </patch_collection>
    <patch_collection type="everywhere">
      <patch ID="1">
        <number>78</number>
        <carveout_patches>
          <patch_collection type="annulus">
            <patch ID="1">
              <x0>0.0</x0>
              <y0>0.0</y0>
              <z0>0.0</z0>
              <inner_radius>100.0</inner_radius>
              <outer_radius>400.0</outer_radius>
            </patch>
          </patch_collection>
          <patch_collection type="elliptical_disc">
            <patch ID="1">
              <x0>-200.0</x0>
              <y0>200.0</y0>
              <z0>0.0</z0>
              <a>100.0</a>
              <b>60.0</b>
              <rotation>pi/6</rotation>
            </patch>
          </patch_collection>
        </carveout_patches>
      </patch>
    </patch_collection>
  </cell_patches>
</ic_cells>
```
