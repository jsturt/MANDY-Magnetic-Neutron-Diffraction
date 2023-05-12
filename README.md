# Magnetic-Neutron-Diffraction (MaND-y)
A python (3.7+) package in development for the simulation of magnetic neutron diffraction from crystals generated with CIF files.

# Requirements
- pandas
- numpy
- matplotlib
- os
- crystals
- gemmi
- requests
- math
- cmath

# Installation
Unless you wish to build the package yourself, pip is the easiest method of installation.
- Navigate to the releases tab on the right-hand side of the page.
- Find the most recent release.
- Under 'Assets' click the `tar.gz` file (i.e. `mandy-X.X.tar.gz`) to begin downloading.
- Navigate to where the file is downloaded and run `pip install`:
```shell
pip install mandy-X.X.tar.gz
```
- You should now be setup to use `mandy`!

# Examples
- Download the `Examples` directory. 
- Contained within are two example projects looking at Spin Density Waves (**SDW**), one for Chromimum (**Cr**) and one of Niobium Iron 2 (**NbFe<sub>2</sub>**)
- Each has .py project files for setting up the structural and magnetic information plus the requirements for the calculations such as:
	- site-names file (Allows association of any magnetic site to any lattice site)
	- CIF files for the description of the unit cells
- The **Cr** example is more simple than the **NbFe<sub>2</sub>** example, so the recommended starting point is `Examples/Cr/Cr_SDW.py`.
