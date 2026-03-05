# Rapidly deployable hulls and on-demand tunable hydrodynamics with shape morphing curved crease origami

## Authors
- Hardik Y. Patil (hardikyp@umich.edu)
- Kevin J. Maki (kjmaki@umich.edu)
- Evgueni T. Filipov (filipov@umich.edu)

This codebase contains:
- the bar-and-hinge formulation for simulating folding of curved-crease origami-inspired hulls;
- automation scripts for the generation of CAD geometry from origami folding simulations; and
- scripts to set up reduced-order hydrodynamic simulations of origami hulls in POWERSEA.

## Abstract
> Traditional hull fabrication relies on labor- and time-intensive methods to generate smooth, curved surfaces. These conventional methods often lead to hull surface topologies that are static in design with hydrodynamics aimed at handling a broad range of sea conditions but not optimized for any specific scenario. In this paper, we introduce a method of rapidly fabricating planing hulls using the principles of curved-crease origami. Starting from a flat-folded state, the curved-crease origami hulls can be deployed to match traditional planing hull shapes like the VPS (deep-V, Planing hull with Straight face) and the GPPH (General Purpose Planing Hull). By extension of the ability to conform to a desired shape, we show that the curved-crease origami hulls can emulate desired hydrodynamic characteristics in still as well as wavy water conditions. Furthermore, we demonstrate the shape-morphing ability of curved-crease origami hulls, enabling them to switch between low and high deadrise configurations. This ability allows for on-demand tuning of the hull hydrodynamic performance. We present proof-of-concept origami hulls to demonstrate the practical feasibility of our method. Hulls fabricated using the curved-crease origami principles can adapt to different sea states, and their flat foldability and deployability facilitate easy transport and deployment for rapid response naval operations such as rescue missions and the launch of crewless aquatic vehicles.

## Citing this work
If you publish results that leverage this codebase or the findings from our paper, please cite the original research article.
> Patil, H. Y., Maki, K. J., & Filipov, E. T. (2024). Rapidly deployable hulls and on-demand tunable hydrodynamics with shape morphing curved crease origami. Journal of Fluids and Structures, 130, 104176. [https://doi.org/10.1016/j.jfluidstructs.2024.104176](https://doi.org/10.1016/j.jfluidstructs.2024.104176)

## PDF
[Read the submitted manuscript here](2024-JLFS-Patil-et-al-open-access).

## Acknowledgement
The bar and hinge formulation in this codebase is modified from Woodruff and Filipov (2020) to simulate origami-inspired ship hull geometries and to perform shape matching by optimizing geometric parameters. Please refer to the paper below for details about the bar and hinge model formulation for curved-crease origami:
> S. R. Woodruff, E. T. Filipov (2020). A bar and hinge model formulation for structural analysis of curved-crease origami. International Journal of Solids and Structures

We would also like to acknowledge the prior work from Ke Liu and Glaucio H. Paulino for straight crease origami, which paved the way for the development of this code.
