## Towards a new standard for seismic moment tensor inversion containing 3D Earth structure uncertainty 

Thanh-Son Pham (thanhson.pham@anu.edu.au)

This repository contains project files to reproduce results and figures reported in "Towards a new standard for seismic moment tensor inversion containing 3D Earth structure uncertainty" by Pham et al. (2024) to be published in Geophysical Journal International.

Content of this package:

- ses3d_r07_b: copy of ses3d program (https://cos.ethz.ch/software/production/ses3d.html) with customized Python script to perturb 3D Earth models and computation of their Green's funcitons.
- MainMTI.py: main python script to excecute primary data analyses reported in the paper.
- Pham_etal_GJI2024.ipynb: Jupyter notebook contains script to reproduce published figures.
- helper.py: Python library of helpful function used in the main scripts.
- IntegrationTests and InversionAssumptionTests contain pre-computed analysis outputs, used to generate the figures.

**Canberra, May 2024.**
