# 3D Printing Airfoil Sections

This repository includes a workflow for turning Eppler/NACA airfoil output into printable STL files.

## How STL files are generated

1. Run the Fortran solver (`bin\profile.exe`) with an input deck from `data\input\`.
   - Example inputs used: `e1098.dat`, `eppler993.dat`, `naca2412.dat`
2. The solver writes `profile.out` in the working directory.
3. Convert `profile.out` to an extruded STL using:

```powershell
python .\src\profile_to_stl.py .\profile.out .\data\output\e1098_section.stl --chord-mm 150 --span-mm 40
```

The converter reads the `N X Y` airfoil coordinate table, scales chord length, extrudes to span, and writes an ASCII STL mesh.

## Where to find generated STL files

Generated airfoil sections are stored in:

- `data\output\e1098_section.stl`
- `data\output\e993_section.stl`
- `data\output\naca2412_section.stl`

## Programs you can use to view STL files

- Paint 3D (Windows)
- Blender (`File -> Import -> STL`)
- MeshLab
- Microsoft 3D Viewer

## How STL is converted to printer-ready format

STL is geometry only; most 3D printers cannot print STL directly.

Use a slicer to convert STL to machine instructions:

1. Open STL in a slicer (PrusaSlicer, Cura, Bambu Studio, OrcaSlicer, etc.).
2. Set units/scale (these files are in mm), orientation, layer height, infill, walls, supports, and material profile.
3. Slice the model to produce printer-ready output (commonly G-code or printer-specific job formats).
4. Send that sliced file to the printer via SD card, USB, LAN, or cloud workflow.

