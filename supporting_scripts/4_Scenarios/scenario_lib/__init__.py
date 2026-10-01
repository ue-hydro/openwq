# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Scenario Library
================

Scenario analysis on a calibrated openWQ model (what-if runs on management,
policy and climate options): the option catalogue and runner (``scenarios``),
the interactive scenario setup report (``Gen_Scenario_Setup_Report``) and the
comparison report (``Gen_Scenario_Results_Report``).

It builds on the calibration framework (``../../3_Calibration/calibration_lib``:
model runner, parameter handler, reach mapping, report helpers), which is put
on ``sys.path`` here when it is not importable yet.
"""
import os as _os
import sys as _sys

try:
    import calibration_lib  # noqa: F401
except ImportError:
    _calib = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "3_Calibration"))
    if _os.path.isdir(_os.path.join(_calib, "calibration_lib")):
        _sys.path.insert(0, _calib)
