# OpenWQ 
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

#### Docker
There is a dockerfile provided in OpenWQ repository. 
 1. cd into the OpenWQ repository (if following the common steps, this would be `summa/build/source/openwq/openwq`)
 2. run `docker build -t openwq .`
 3. cd back up to the summa directory
 4. run `docker run -itd -v $(pwd):/code openwq`

#### Apptainer

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
* Used for post-processing
* See [supporting_scripts](supporting_scripts) folder

## Working Example
* Coming soon!
<!--
* See ["Working_example"](Working_example) folder
* Updates coming soon!
-->
