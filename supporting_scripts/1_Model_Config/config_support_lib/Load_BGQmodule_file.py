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
# !/usr/bin/env python3
"""
OpenWQ JSON File Copy with Header Comments
Copy JSON file and add header comments if they don't already exist
"""

import json
import re
from pathlib import Path
from typing import Dict, List, Optional


def _apply_global_parameter_overrides(source_content: str, overrides: Dict[str, float]) -> str:
    """Set values of the template's GLOBAL_PARAMETERS block in the copied text.

    Only the numbers are replaced, so the template keeps its layout. Every key must
    already exist in GLOBAL_PARAMETERS (a typo would otherwise be silently ignored).
    """
    try:
        gp = json.loads(source_content).get("GLOBAL_PARAMETERS")
    except json.JSONDecodeError as e:
        raise ValueError(f"bgc_global_parameter_overrides: the BGC template is not plain JSON ({e})")
    if not isinstance(gp, dict):
        raise ValueError("bgc_global_parameter_overrides: the BGC template has no GLOBAL_PARAMETERS block")
    unknown = [k for k in overrides if k not in gp]
    if unknown:
        raise ValueError(f"bgc_global_parameter_overrides: {unknown} not in the template's GLOBAL_PARAMETERS")
    i0 = source_content.index('"GLOBAL_PARAMETERS"')
    i1 = source_content.index("}", i0)          # the block holds numbers only, no nested objects
    block = source_content[i0:i1]
    for key, value in overrides.items():
        block, n = re.subn(r'("%s"\s*:\s*)[-+0-9.eE]+' % re.escape(key),
                           lambda m: m.group(1) + repr(float(value)), block)
        if n != 1:
            raise ValueError(f"bgc_global_parameter_overrides: '{key}' found {n} times in GLOBAL_PARAMETERS")
    out = source_content[:i0] + block + source_content[i1:]
    new = json.loads(out)["GLOBAL_PARAMETERS"]
    for key, value in overrides.items():
        print(f"  GLOBAL_PARAMETERS override: {key} = {new[key]} (template {gp[key]})")
    return out


def load_bgq_module_json(
        json_header_comment: List[str],
        path2selected_NATIVE_BGC_FLEX_framework: str,
        bgc_config_filepath: str,
        global_parameter_overrides: Optional[Dict[str, float]] = None
) -> None:
    """
    Copy a JSON file from source to destination and add header comments if they don't exist.

    Parameters:
        json_header_comment: List of comment lines to add at the top of the file
        path2selected_NATIVE_BGC_FLEX_framework: Source file path to copy from
        bgc_config_filepath: Destination file path to copy to
        global_parameter_overrides: optional {name: value} written into the copy's
            GLOBAL_PARAMETERS, for site-specific constants of a shared template
            (e.g. the water-temperature relation a_Tw, b_Tw of each basin)
    """

    # Create the destination directory path if it doesn't exist
    dest_path = Path(bgc_config_filepath)
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    # Read the source file
    source_path = Path(path2selected_NATIVE_BGC_FLEX_framework)

    if not source_path.exists():
        raise FileNotFoundError(f"Source file not found: {path2selected_NATIVE_BGC_FLEX_framework}")

    with open(source_path, 'r') as f:
        source_content = f.read()

    # Strip any leading/trailing whitespace from source content
    source_content = source_content.strip()

    if global_parameter_overrides:
        source_content = _apply_global_parameter_overrides(source_content, global_parameter_overrides)

    # Check if destination file already exists and has comments
    has_comments = False
    if dest_path.exists():
        with open(bgc_config_filepath, 'r') as f:
            dest_content = f.read()
        # Check if any of the header comments already exist in the destination file
        for comment in json_header_comment:
            if comment.strip() in dest_content:
                has_comments = True
                break

    # Write to destination file
    with open(bgc_config_filepath, 'w') as f:
        # Add header comments if they don't exist and are provided
        if json_header_comment and not has_comments:
            for comment in json_header_comment:
                f.write(comment + "\n")
        # Write the source content
        f.write(source_content)
        # Add newline at end of file
        f.write("\n")

    print(f"  Biogeochemistry cycling framework file saved to: {bgc_config_filepath}")
