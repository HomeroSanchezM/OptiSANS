SANS Deuteration Optimization
==============================

A genetic algorithm framework for finding the optimal deuteration pattern of a
protein for Small-Angle Neutron Scattering (SANS) contrast variation
experiments.

The algorithm evolves a population of chromosomes, each representing which
amino acid types are fully deuterated and the D2O percentage of the solvent.
Each chromosome is used to generate a deuterated PDB file, which is then
passed to a SANS simulation program (Pepsi-SANS). The fitness of a solution
is defined as the product of the areas between the scaled simulated curve and
two experimental references, multiplied by the signal-to-background ratio.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   configuration
   visualization
   api
   
