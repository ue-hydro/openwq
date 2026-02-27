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
Reach Mapping Module
====================

Maps external reach IDs (reachID, hruId) to internal OpenWQ indices (ix, iy, iz).

This allows users to specify observation locations and source/sink loads using
the familiar reach IDs from their model rather than internal array indices.
"""

import h5py
import numpy as np
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import logging

logger = logging.getLogger(__name__)


class ReachMapper:
    """
    Maps between external reach IDs and internal OpenWQ indices.

    The mapping is extracted from HDF5 output files where OpenWQ stores:
    - 'reachID' (or 'hruId'): External identifiers
    - 'xyz_elements': Internal (ix, iy, iz) indices

    This allows users to specify locations using reach IDs instead of
    having to know the internal array indices.
    """

    def __init__(self,
                 mapping_source: Union[str, Path] = None,
                 mapping_key: str = "reachID"):
        """
        Initialize reach mapper.

        Parameters
        ----------
        mapping_source : str or Path
            Path to either:
            - An HDF5 output file containing the mapping
            - A directory containing HDF5 files (will use first valid file)
            - A JSON mapping file
        mapping_key : str
            The mapping key name ("reachID" for mizuRoute, "hruId" for SUMMA)
        """
        self.mapping_key = mapping_key
        self.reach_to_xyz: Dict[str, Tuple[int, int, int]] = {}
        self.xyz_to_reach: Dict[Tuple[int, int, int], str] = {}

        if mapping_source:
            self.load_mapping(mapping_source)

    def load_mapping(self, source: Union[str, Path]) -> bool:
        """
        Load the reach ID to xyz mapping from source.

        Parameters
        ----------
        source : str or Path
            Path to mapping source (HDF5 file, directory, or JSON)

        Returns
        -------
        bool
            True if mapping was loaded successfully
        """
        source = Path(source)

        if source.is_file():
            if source.suffix in ['.h5', '.hdf5']:
                return self._load_from_hdf5(source)
            elif source.suffix == '.json':
                return self._load_from_json(source)
        elif source.is_dir():
            return self._load_from_directory(source)

        logger.warning(f"Could not load mapping from: {source}")
        return False

    def _load_from_hdf5(self, h5_file: Path) -> bool:
        """Load mapping from an HDF5 file."""
        try:
            with h5py.File(h5_file, 'r') as f:
                # Find the mapping key
                if self.mapping_key in f:
                    reach_ids = f[self.mapping_key][:]
                elif 'reachID' in f:
                    reach_ids = f['reachID'][:]
                    self.mapping_key = 'reachID'
                elif 'hruId' in f:
                    reach_ids = f['hruId'][:]
                    self.mapping_key = 'hruId'
                else:
                    logger.warning(f"No mapping key found in {h5_file}")
                    return False

                # Get xyz_elements
                if 'xyz_elements' not in f:
                    logger.warning(f"No xyz_elements found in {h5_file}")
                    return False

                xyz = f['xyz_elements'][:]

                # Build mapping
                # xyz shape is (3, n_elements) where rows are [ix, iy, iz]
                for i in range(len(reach_ids)):
                    # Decode bytes to string if needed
                    reach_id = reach_ids[i]
                    if isinstance(reach_id, bytes):
                        reach_id = reach_id.decode('utf-8')
                    reach_id = str(reach_id)

                    # Get indices (convert from 1-based to 0-based if needed)
                    ix = int(xyz[0, i])
                    iy = int(xyz[1, i])
                    iz = int(xyz[2, i])

                    self.reach_to_xyz[reach_id] = (ix, iy, iz)
                    self.xyz_to_reach[(ix, iy, iz)] = reach_id

                logger.info(f"Loaded mapping for {len(self.reach_to_xyz)} elements from {h5_file}")
                return True

        except Exception as e:
            logger.error(f"Error loading HDF5 mapping: {e}")
            return False

    def _load_from_directory(self, directory: Path) -> bool:
        """Load mapping from first valid HDF5 file in directory."""
        # Look for HDF5 files
        h5_patterns = ['*.h5', '*.hdf5']

        for pattern in h5_patterns:
            for h5_file in directory.glob(pattern):
                if self._load_from_hdf5(h5_file):
                    return True

        # Try subdirectories (e.g., HDF5/)
        for subdir in directory.iterdir():
            if subdir.is_dir():
                for pattern in h5_patterns:
                    for h5_file in subdir.glob(pattern):
                        if self._load_from_hdf5(h5_file):
                            return True

        logger.warning(f"No valid HDF5 files found in {directory}")
        return False

    def _load_from_json(self, json_file: Path) -> bool:
        """Load mapping from a JSON file."""
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)

            # Expected format: {"reach_id": [ix, iy, iz], ...}
            for reach_id, xyz in data.items():
                reach_id = str(reach_id)
                ix, iy, iz = int(xyz[0]), int(xyz[1]), int(xyz[2])
                self.reach_to_xyz[reach_id] = (ix, iy, iz)
                self.xyz_to_reach[(ix, iy, iz)] = reach_id

            logger.info(f"Loaded mapping for {len(self.reach_to_xyz)} elements from {json_file}")
            return True

        except Exception as e:
            logger.error(f"Error loading JSON mapping: {e}")
            return False

    def save_mapping(self, json_file: Path):
        """
        Save the current mapping to a JSON file.

        Useful for creating a static mapping file that can be reused.

        Parameters
        ----------
        json_file : Path
            Output JSON file path
        """
        # Convert tuple keys to list for JSON serialization
        data = {reach_id: list(xyz) for reach_id, xyz in self.reach_to_xyz.items()}

        with open(json_file, 'w') as f:
            json.dump(data, f, indent=2)

        logger.info(f"Saved mapping to {json_file}")

    def get_xyz(self, reach_id: Union[str, int]) -> Optional[Tuple[int, int, int]]:
        """
        Get internal (ix, iy, iz) indices for a reach ID.

        Parameters
        ----------
        reach_id : str or int
            External reach ID

        Returns
        -------
        Tuple[int, int, int] or None
            (ix, iy, iz) tuple, or None if not found
        """
        reach_id = str(reach_id)
        return self.reach_to_xyz.get(reach_id)

    def get_reach_id(self, ix: int, iy: int, iz: int) -> Optional[str]:
        """
        Get reach ID for internal indices.

        Parameters
        ----------
        ix, iy, iz : int
            Internal indices

        Returns
        -------
        str or None
            Reach ID, or None if not found
        """
        return self.xyz_to_reach.get((ix, iy, iz))

    def get_ix(self, reach_id: Union[str, int]) -> Optional[int]:
        """Get ix index for a reach ID."""
        xyz = self.get_xyz(reach_id)
        return xyz[0] if xyz else None

    def convert_reach_ids_to_xyz(self,
                                  reach_ids: List[Union[str, int]]
                                  ) -> List[Tuple[int, int, int]]:
        """
        Convert a list of reach IDs to xyz indices.

        Parameters
        ----------
        reach_ids : List
            List of reach IDs

        Returns
        -------
        List[Tuple[int, int, int]]
            List of (ix, iy, iz) tuples
        """
        result = []
        for rid in reach_ids:
            xyz = self.get_xyz(rid)
            if xyz:
                result.append(xyz)
            else:
                logger.warning(f"Reach ID not found in mapping: {rid}")
        return result

    def is_valid_reach_id(self, reach_id: Union[str, int]) -> bool:
        """Check if a reach ID exists in the mapping."""
        return str(reach_id) in self.reach_to_xyz

    def get_all_reach_ids(self) -> List[str]:
        """Get list of all mapped reach IDs."""
        return list(self.reach_to_xyz.keys())

    def get_mapping_summary(self) -> Dict:
        """Get summary of the mapping."""
        if not self.reach_to_xyz:
            return {"status": "empty", "n_elements": 0}

        # Get range of indices
        all_ix = [xyz[0] for xyz in self.reach_to_xyz.values()]
        all_iy = [xyz[1] for xyz in self.reach_to_xyz.values()]
        all_iz = [xyz[2] for xyz in self.reach_to_xyz.values()]

        return {
            "status": "loaded",
            "n_elements": len(self.reach_to_xyz),
            "mapping_key": self.mapping_key,
            "ix_range": (min(all_ix), max(all_ix)),
            "iy_range": (min(all_iy), max(all_iy)),
            "iz_range": (min(all_iz), max(all_iz)),
            "sample_mappings": dict(list(self.reach_to_xyz.items())[:5])
        }

    def __len__(self):
        """Return number of mapped elements."""
        return len(self.reach_to_xyz)

    def __contains__(self, reach_id):
        """Check if reach ID is in mapping."""
        return str(reach_id) in self.reach_to_xyz

    def __repr__(self):
        return f"ReachMapper(n_elements={len(self)}, mapping_key='{self.mapping_key}')"


def create_mapping_from_output(output_dir: Union[str, Path],
                               output_json: Union[str, Path] = None,
                               mapping_key: str = "reachID") -> ReachMapper:
    """
    Convenience function to create a ReachMapper from OpenWQ output directory.

    Parameters
    ----------
    output_dir : str or Path
        Path to openwq_out directory
    output_json : str or Path, optional
        If provided, save the mapping to this JSON file for reuse
    mapping_key : str
        Mapping key name ("reachID" or "hruId")

    Returns
    -------
    ReachMapper
        Initialized mapper with loaded data

    Example
    -------
    >>> mapper = create_mapping_from_output("/path/to/openwq_out")
    >>> xyz = mapper.get_xyz("1200014181")
    >>> print(f"Reach 1200014181 is at index {xyz}")
    """
    mapper = ReachMapper(mapping_key=mapping_key)
    mapper.load_mapping(Path(output_dir))

    if output_json and len(mapper) > 0:
        mapper.save_mapping(Path(output_json))

    return mapper


def convert_ss_json_reach_to_xyz(ss_json_path: Union[str, Path],
                                  mapper: ReachMapper,
                                  output_path: Union[str, Path] = None) -> Dict:
    """
    Convert reach IDs in a Source/Sink JSON file to xyz indices.

    This allows users to write SS files using reach IDs, then convert
    them to the internal format before running the model.

    Parameters
    ----------
    ss_json_path : str or Path
        Path to source/sink JSON file
    mapper : ReachMapper
        Initialized reach mapper
    output_path : str or Path, optional
        If provided, write converted JSON to this path

    Returns
    -------
    Dict
        Converted JSON data

    Example
    -------
    >>> mapper = create_mapping_from_output("/path/to/openwq_out")
    >>> convert_ss_json_reach_to_xyz(
    ...     "openWQ_fertilizer_N.json",
    ...     mapper,
    ...     "openWQ_fertilizer_N_converted.json"
    ... )
    """
    with open(ss_json_path, 'r') as f:
        # Read with comment handling
        content = f.read()
        # Remove // comments (simple approach)
        lines = []
        for line in content.split('\n'):
            if '//' in line:
                line = line[:line.index('//')]
            lines.append(line)
        content = '\n'.join(lines)
        data = json.loads(content)

    # Process each entry
    for key, entry in data.items():
        if key == "METADATA":
            continue

        if "Data" not in entry:
            continue

        # Process data entries
        new_data = {}
        for data_key, data_val in entry["Data"].items():
            if not isinstance(data_val, list):
                new_data[data_key] = data_val
                continue

            # Check if ix, iy, iz are reach IDs
            # Format: [YYYY, MM, DD, HH, min, sec, ix, iy, iz, load, type]
            if len(data_val) >= 10:
                ix_val = data_val[6]
                iy_val = data_val[7]
                iz_val = data_val[8]

                # If ix is a string (reach_id), convert it
                if isinstance(ix_val, str) and ix_val != "all":
                    xyz = mapper.get_xyz(ix_val)
                    if xyz:
                        new_val = data_val.copy()
                        new_val[6] = xyz[0]
                        new_val[7] = xyz[1]
                        new_val[8] = xyz[2]
                        new_data[data_key] = new_val
                        logger.info(f"Converted reach {ix_val} to ix={xyz[0]}, iy={xyz[1]}, iz={xyz[2]}")
                    else:
                        logger.warning(f"Could not find mapping for reach ID: {ix_val}")
                        new_data[data_key] = data_val
                else:
                    new_data[data_key] = data_val
            else:
                new_data[data_key] = data_val

        entry["Data"] = new_data

    # Write output if path provided
    if output_path:
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        logger.info(f"Wrote converted SS file to {output_path}")

    return data
