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

from pathlib import Path
from typing import List


def load_bgq_module_json(
        json_header_comment: List[str],
        path2selected_NATIVE_BGC_FLEX_framework: str,
        bgc_config_filepath: str
) -> None:
    """
    Copy a JSON file from source to destination and add header comments if they don't exist.

    Parameters:
        json_header_comment: List of comment lines to add at the top of the file
        path2selected_NATIVE_BGC_FLEX_framework: Source file path to copy from
        bgc_config_filepath: Destination file path to copy to
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

    print(f"✓ Biogeochemistry cycling framework file saved to: {bgc_config_filepath}")