# AutoCAD Script Documentation

(*In development. Last updated July 8, 2022*)

This document describes the procedure followed to export nodal coordinates obtained from the folded shape using Bar and Hinge model and generate a 3D geometry using AutoCAD for hydrodynamic simulations. The nodal coordinates are organised in a structure (MATLAB datatype) called ``curves``. The curve name and its corresponding index in the structure is shown below

## ``curves``

1. Starboard side chine
2. Keel
3. Port side chine
4. Port side top edge
5. Starboad side top edge
6. Stern edge between keel and starboard side chine
7. Stern edge between keel and port side chine
8. Stern edge between starboard side chine and top edge
9. Stern edge between port side chine and top edge
10. Stern edge between startboard and port side top edges

## Procedure

1. Draw all curves using the ``SPLINE`` function in AutoCAD
2. Loft surfaces between adjacent splines using the ``LOFT`` function
3. Convert the enclosed volume into a 3D shape using ``SURFSCULPT`` function
4. Rotate and align the 3D solid such that the boat faces negative X direction
5. Use the function ``MASSPROP`` to find geometric properties of the boat like center of gravity, radii of gyration, etc.
6. Delete the 3D solid and keep only the *Keel* and *Starboard Chine* curves
7. Export the curves in an ``.igs`` file format
