# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
OpenWQ Calibration Framework
============================

A comprehensive calibration framework for the OpenWQ water quality model.

Features:
- DDS (Dynamically Dimensioned Search) optimization
- Morris and Sobol sensitivity analysis
- Support for Docker and Apptainer (HPC) execution
- Checkpointing and restart capability
- Multi-objective calibration support

Usage:
    1. Copy calibration_config_template.py to your working directory
    2. Edit the configuration with your parameters and bounds
    3. Run: python your_calibration_config.py

Documentation: https://openwq.readthedocs.io
"""

__version__ = "1.0.0"
__author__ = "Diogo Costa"

from .calibration_lib.calibration_driver import run_calibration, run_sensitivity_analysis
