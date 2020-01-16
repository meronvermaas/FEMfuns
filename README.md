<p align="center">
  <img src="https://raw.githubusercontent.com/meronvermaas/FEMfuns/master/logo.png" width="250">
</p>

Welcome to FEMfuns (Finite Element Method for useful neuroscience simulations).

FEMfuns allows neuroscientists to solve the forward problem in a variety of different geometrical domains, including various types of source models and electrode properties, such as resistive and capacitive materials.

The results of using FEMfuns demonstrate that the pipeline presented in this paper is an accurate and fexible tool to simulate signals generated on electrode grids by the spatiotemporal electrical activity patterns produced by cortical neurons and thereby allows the user to optimize grids for brain computer interfaces including exploration of alternative electrode materials/properties.

This software is provided under the GNU General Purpose License version 3.0, You will find a copy of this license within this folder, or from online here: https://www.gnu.org/licenses/gpl-3.0.txt

Instructions to run code:

install anaconda https://docs.continuum.io/anaconda/install/

clone github code https://github.com/meronvermaas/FEMfuns

move to cloned FEMfuns directory: cd FEMfuns

conda env create -f environment.yml
conda activate femfuns
conda develop pipeline_code/

Get the mesh geometries

Run the code and wait (depending on the study and geometry) a while:

Study 1 run the following lines in terminal:
python3.7 fem_study1.py
python3.7 analytical_correct.py
python3.7 plot_study1.py

Study 2 run the following lines in terminal:
python3.7 fem_study2.py
python3.7 plot_study2.py

Study 3 run the following lines in terminal:
python3.7 fem_study3.py

Inspect pvd files from study3 in Paraview (https://www.paraview.org/download/)
