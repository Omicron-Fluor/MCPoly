# MCPoly
Some methods to deal with some manipulation of computational chemistry, mostly about mechanical properties of polymers.
![](./reference/function.png)

## Overview
`MCPoly` is a Python library to make some steps of computational chemistry easier. It includes some functions of drawing molecule structures, creating proper .xyz , .inp , .mol and .data file, which specialised for using ORCA and LAMMPS, and some functions for researching the mechanical property of some polymers.

## Updated in v0.4.1 (29.05.23)
1. BIG UPDATE: Molecule Designer GUI is included in this package. You can now use it to draw molecule with 3D view at the same time.
<img src="https://github.com/Omicron-Fluor/MCPoly/blob/main/reference/moldraw_gui.png">
2. Add 'chain' in 'lmpset' package. You can now use it to convert periodic polymer structure input data now.

## Updated in v0.3.0 (19.05.23)
1. BIG UPDATE: Now status package can be used to see energy and gibbs free energy by .energy and .gibbs.
2. BIG UPDATE: The new Reaction Diagram GUI is included in this package. You can use it to merge some tables and sea the relevant reaction diagram as a reference.
<img src="https://github.com/Omicron-Fluor/MCPoly/blob/main/reference/status_gui.png">
3. Fix some bugs about availability of Stress-Strain Curve Output.

## Updated in v0.2.1 (16.05.23)
1. BIG UPDATE: Add GUI ORCA setting screens for multiple ORCA calculations and Stress-Strain Curve Output. You can use orca.mgui and sscurve.gui to learn more.
2. BIG UPDATE: MCPoly.vis is created. Use MCPoly.vis() to directly show all four GUI screen.
<img src="https://github.com/Omicron-Fluor/MCPoly/blob/main/reference/vis.png">
3. Fix some bugs on mechanical calculations and stress-strain curve drawings.

## Functions for ORCA
<img src="https://github.com/Omicron-Fluor/MCPoly/blob/main/reference/ORCA.png" width="400" height="263" >

### orcaset
Used to create ORCA input files and run it on ORCA. It's especially handy for researching mechanical property of polymers.

### status
Because we can't use ORCA to visualize the geometry structure, this command can be used to see the optimisation status and the trait of geometry structure.

### view3d
See the 3D structure of a normal .xyz file.

### sscurve
With calculated .xyz file, we can draw the stress-strain curve of each polymer, and we can also calculate the Young's modulus of relevant polymers.

### moldraws (Under Construction)
Used to build a simple molecule and save it under .xyz form.

## Functions for LAMMPS
<img src="https://github.com/Omicron-Fluor/MCPoly/blob/main/reference/LAMMPS.png" width="395" height="100" >

### lmpset
Used to draw special patterns of polymers. Mostly in grids.

## Installation
To get `MCPoly`, you can install it with pip:
    `$ pip install MCPoly`

If you want to get the latest version of `MCPoly`, you can see the latest release here:

<https://github.com/Omicron-Fluor/MCPoly/release> 

There will be a corresponding release on `pip` for each release on GitHub, and you can update your `MCPoly` with:

`$ pip install MCPoly --upgrade`

## How to cite
<https://github.com/Omicron-Fluor/MCPoly>

## Outlook
We will add some new function about polymers based on ReaxFF.
