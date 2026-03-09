# Guide

You can create a simple XML by running
```julia
using PhysiCellCellCreator
path_to_folder = "template_folder"
createICCellXMLTemplate(path_to_folder) # will create cells.xml inside the folder "template_folder"
```

Then, you can edit the XML file as needed.
Note, the root of the XML file is `<ic_ecm>`.

## Cell patches
The children of the root are `<cell_patches>` elements, each corresponding to a different cell type in the PhysiCell simulation.
The cell type is set by the attribute `name`, e.g., `<cell_patches name="fibroblast">`.

## Patch collection
Each `<cell_patches>` element contains a collection of `<patch_collection>` elements, each containing a particular type of patch.
The currently supported patch types are `everywhere`, `disc`, `annulus`, and `rectangle`.
A typical example is `<patch_collection type="everywhere">`.

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

A `disc` patch can also have a `normal` element, defining the normal vector of the disc.
If omitted, the normal vector is assumed to be `(0, 0, 1)`.
To set this, add a child element to the patch as `<normal>1.0,0.0,0.0</normal>`.
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

## Refined control of cell placement

### Restrict to domain
Each of the above patches (except `everywhere`) can have a `restrict_to_domain` element, which will force cells to be placed within the domain as specified by the `domain_dict` parameter in the `generateICCell` function.
If omitted, it is assumed to be "true".
Possible values are "true", "false", "0", and "1".
Add the following to the patch element to set: `<restrict_to_domain>false</restrict_to_domain>`.

### Carveouts
Every patch type can also include a `<carveout_patches>` element, which specifies regions where cells should not be placed.
The currently supported carveouts are `disc`, `annulus`, and `rectangle`.
These follow the same structure as the patches above, but do not contain a `number` element.
`disc` and `annulus` carveouts cannot have a `normal` element, they are assumed to be parallel to the xy-plane.

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
          <patch_collections type="disc">
              <patch ID="1">
                  <x0>130.0</x0>
                  <y0>130.0</y0>
                  <z0>0.0</z0>
                  <radius>15.0</radius>
              </patch>
          </patch_collections>
      </carveout_patches>
  </patch>
</patch_collection>
```

### Max fails
When attempting to place cells, points taken uniformly randomly from the entire patch are selected and then pared down based on the domain and careveouts.
This can lead to a situation where many samplings must be taken to finally achieve `number` of cells in each patch.
In the worst case scenario, a patch has no area and cannot support any cells.
To prevent the program from running indefinitely, a `max_fails` parameter controls how many samples are taken before giving up.
If omitted, it is assumed to be 100.
Add the following to the patch element (not the carveout, but the patch itself) to set: `<max_fails>100</max_fails>`.

## Example XML
A full example XML is shown below.
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
        <width>100.0</width>
        <height>100.0</height>
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
        </carveout_patches>
      </patch>
    </patch_collection>
  </cell_patches>
</ic_cells>
```