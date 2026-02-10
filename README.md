<p align="left">
  <img src="uevora_logo.png" alt="University of Évora" height="80">
</p>

# OpenWQ

**A flexible biogeochemical modeling framework for water quality simulations**

OpenWQ couples with existing hydrological and hydrodynamic models to enable multi-scale, multi-chemistry simulations. It supports both kinetic reaction networks (via the native BGC-Flex engine) and equilibrium geochemistry (via PHREEQC integration).

[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://openwq.readthedocs.io)
[![License](https://img.shields.io/badge/license-GPL--3.0-green)](LICENSE)

---

## Table of Contents
- [OpenWQ](#openwq)
	- [Table of Contents](#table-of-contents)
	- [Introduction](#introduction)
	- [Branches](#branches)
		- [Active](#active)
		- [Notes](#notes)
	- [Compiling](#compiling)
		- [Summa-OpenWQ](#summa-openwq)
			- [Common Steps](#common-steps)
			- [macOS](#macos)
			- [Linux](#linux)
			- [Docker](#docker)
			- [Apptainer](#apptainer)
	- [Execution](#execution)
	- [Visualization](#visualization)
	- [Supporting Scripts](#supporting-scripts)
	- [Working Example](#working-example)

## Introduction
* Soure code for the multi-chemistry, multi-scale OpenWQ modelling framework.
* Reading material:
	* [Documentation](https://openwq.readthedocs.io)
	* Articles: Coming soon!
<!-- Sample syntax from FLUXOS README
	* Theoretical background (original FLUXOS):
		* [EMS paper](https://www.sciencedirect.com/science/article/pii/S1364815216306193?via%3Dihub)
		* [PhD thesis](https://scholarbank.nus.edu.sg/handle/10635/124183)
	* Applications (original FLUXOS):
		* [STC paper](https://www.sciencedirect.com/science/article/pii/S0169772216300948?via%3Dihub)
		* [JCH paper](https://www.sciencedirect.com/science/article/pii/S0169772216300948?via%3Dihub)
		* [JAWRA](https://onlinelibrary.wiley.com/doi/full/10.1111/1752-1688.12316)
	* FLUXOS-OVERLAND
		* [Poster](https://www.researchgate.net/publication/333324452_Hydrodynamic_modelling_of_snowmelt_flooding_events_and_nutrient_transport_in_the_Canadian_Prairies_using_the_FLUXOS_model?channel=doi&linkId=5ce70f0a458515712ebda98b&showFulltext=true)
-->

## Branches
### Active
* main: primary branch with latest verified updates to the working code 
* development: used to verify updates from feature branches before merging with main
* supporting_scripts: feature branch for development of supporting Python and MATLAB scripts
* wrapper_interfaces: feature branch for development of wrapper interfaces between SUMMA and OpenWQ, and MESH and OpenWQ
* expression_evaluation_optimization: feature branch for optimizing string parsing and the computation of mathematical expressions
* optimization_sink_source: feature branch

### Notes
* Note that the primary branch (formerly known as "master") was renamed as "main"
* If using a local clone with the previous naming scheme, the clone may be updated using the following Git commands:
```
$ git branch -m master main
$ git fetch -p origin
$ git branch -u origin/main main
```

## Compiling
* CMake: [CMakeLists.txt](CMakeLists.txt) is provided
* Library dependencies: Armadillo, VTK, OpenMP 
* CMake minimum version: 3.10

### Summa-OpenWQ
Below are the steps to compile SUMMA with OpenWQ integration.

#### Common Steps
 1. Clone SUMMA: `git clone -b develop https://github.com/ashleymedin/summa.git`
 2. `cd summa/build/source/openwq`
 3. Clone OpenWQ: `git clone -b develop git@github.com:ue-hydro/openwq.git`
   
#### macOS
Native macOS compilation is currently not supported. Please use Docker.

#### Linux
  1. Install or Load (Cluster) Dependencies
    - NetCDF
    - HDF5
    - OpenBLAS or LAPACK
    - Armadillo
      - Included install script in `summa/build/source/openwq/openwq/utils/install_armadillo.sh`
      - Running the script will install Armadillo in the `summa/build/source/openwq/openwq/utils/armadillo-VERSION` directory
    - C++ Compiler (e.g., g++)
    - Fortran Compiler (e.g., gfortran)
    - CMake
  2. Method 1: `cd summa/build/source/openwq/openwq/build/`
  3. `cmake ..`
  4. `make -j 2`
  5. Method 2: `cd summa/build/cmake/`
  6. Edit the `build.pc.cmake` file to include the OpenWQ library with `-DUSE_OPENWQ=ON`
  7. `./build.pc.cmake`

#### Docker
Container files (Dockerfile, docker-compose.yml) are located in the `containers/` folder.

 1. cd into the OpenWQ repository (if following the common steps, this would be `summa/build/source/openwq/openwq`)
 2. cd into containers: `cd containers`
 3. Build using docker-compose: `docker-compose build`
 4. Start the container: `docker-compose up -d`
 5. Connect to the container with `docker attach docker_openwq` or VSCode's Remote-Containers extension
 6. cd into `/code/build/cmake` and adjust the `build.pc.cmake` file to include the OpenWQ library with `-DUSE_OPENWQ=ON`
 7. run `./build.pc.cmake` to compile SUMMA with OpenWQ
 8. Alternatively, you can `cd /code/build/source/openwq/openwq`
 9. Edit the `CMakeLists.txt`:
    - Ensure that `COMPILE_TARGET` is set to `summa_openwq`
    - To use sundials ensure `SOLVER_TYPE` is set to `sundials`, otherwise use `None`
 10. mkdir build && cd build
 11. cmake ..
 12. make -j 2

Alternative (without docker-compose):
 1. `cd containers`
 2. `docker build -t openwq .`
 3. cd back up to the parent directory (e.g., summa)
 4. `docker run -itd -v $(pwd):/code openwq`
 5. Connect to the container with `docker attach <container_id>`

#### Apptainer
The Apptainer definition file `openwq_apptainer.def` is located in the `containers/` folder. This file can be used to build an Apptainer (Singularity) container. NOTE: You will need sudo access to build the container. Once the container is built you can transfer the container to a system with Apptainer to run the container. Running the container does not require sudo access.

TO BUILD THE CONTAINER:
 1. Follow the common steps above to obtain the source code
 2. `cd summa/build/source/openwq/openwq/containers`
 3. `sudo apptainer build openwq.sif openwq_apptainer.def`

COMPILE SUMMA-OPENWQ WITH THE CONTAINER:
 1. There is a script in `utils/` called `compile_summa_apptainer.sh` that will compile SUMMA-OpenWQ using the container.
 2. `cd utils`
 3. `./compile_summa_apptainer.sh`
		- This script will launch the container and compile SUMMA-OpenWQ
		- This will create a summa-openwq that can then be run within the container.

RUNNING SUMMA-OPENWQ WITH THE CONTAINER:
 1. There is an example script in `utils/` called `run_summa_apptainer.sh` that will run SUMMA-OpenWQ using the container.
 2. `cd utils`
 3. `./run_summa_apptainer.sh`
	  - This script will launch the container and run the summa-openwq executable
	  - This will run the summa-openwq executable within the container.
	  - This script depends on the synthetic_tests: https://github.com/KyleKlenk/synthetic_tests, which will need to be cloned into `summa`
 

## Execution
* Coming soon!
<!--
* Create a folder with name "Results" inside the working directory where the input files and executable are located
* input files (see example in Working_example folder)
	* primary input file: e.g., modset
	* DEM file (Esri ASCII-format raster with headers removed ->  this will be fixed soon)
	* DEM of the basin (sub-set of the main DEM file for FLUXOS to know where the boundaries of the basin are)
	* Snowmelt timeseries (time,mm/day)
* to execute: ./fluxos_cpp "argument_1" (where "argument_1" is the primary input file)
-->

## Visualization
* Coming soon!
<!--
* Output stored in "Results" folder may be visualized using [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/):  
[![alt text](https://wci.llnl.gov/sites/wci/files/visit-home.jpg "VisIt")](https://wci.llnl.gov/simulation/computer-codes/visit/)
-->

## Supporting Scripts
* Used for post-processing, calibration, and model configuration
* See [supporting_scripts](supporting_scripts) folder

### Python Environment Setup

The supporting scripts require Python 3.8+ with several dependencies. A virtual environment setup is provided:

```bash
cd supporting_scripts

# Automated setup (recommended)
./setup_venv.sh

# Activate
source .venv/bin/activate

# Or use conda
conda create -n openwq python=3.10
conda activate openwq
pip install -r requirements.txt
```

### Calibration and Sensitivity Analysis
OpenWQ includes a comprehensive calibration framework with:
- **106+ Calibratable Parameters** across all modules (BGC, sorption, sediment, transport, etc.)
- **Sensitivity Analysis**: Morris screening and Sobol variance-based methods
- **DDS Optimization**: Efficient parameter calibration for expensive models
- **HPC Support**: SLURM/PBS templates for high-performance computing

Key documentation:
- [Sensitivity Analysis Guide](supporting_scripts/Calibration/SENSITIVITY_ANALYSIS_GUIDE.md) - Complete parameter sensitivity reference
- [Parameter Defaults](supporting_scripts/Calibration/calibration_lib/parameter_defaults.py) - Scientifically-grounded parameter ranges

## PHREEQC Geochemistry Module

OpenWQ supports the PHREEQC geochemical engine for advanced water quality simulations including:
- Aqueous speciation and saturation indices
- Equilibrium with mineral phases and gases
- Ion exchange and surface complexation
- Kinetic reactions (dissolution, precipitation, redox)

For setup instructions, see the [Documentation](https://openwq.readthedocs.io).

## Working Example
* Coming soon!
<!--
* See ["Working_example"](Working_example) folder
* Updates coming soon!
-->
