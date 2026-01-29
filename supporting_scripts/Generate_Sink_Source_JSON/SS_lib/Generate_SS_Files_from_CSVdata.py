# Copyright 2020, Diogo Costa, diogo.costa@uevora.pt
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

import json
import numpy as np
from pathlib import Path
from typing import List, Dict, Any


def compact_json_dump(data: Dict[str, Any], indent: int = 4) -> str:
    """
    Format JSON with arrays on single lines for better readability.

    Args:
        data: Dictionary to format as JSON
        indent: Number of spaces for indentation

    Returns:
        Formatted JSON string with compact arrays
    """

    def format_value(obj: Any, level: int = 0) -> str:
        indent_str = ' ' * (indent * level)

        if isinstance(obj, dict):
            if not obj:
                return '{}'
            items = [
                f'{indent_str}{json.dumps(k)}: {format_value(v, level + 1)}'
                for k, v in obj.items()
            ]
            return '{\n' + ',\n'.join(items) + '\n' + ' ' * (indent * (level - 1)) + '}'

        elif isinstance(obj, list):
            # Keep arrays on one line
            return '[' + ', '.join(json.dumps(item) for item in obj) + ']'

        else:
            return json.dumps(obj)

    return format_value(data, 1)


def validate_input_consistency(
        chemNames: List[str],
        comptNames: List[str],
        ssType: List[str],
        ssUnits: List[str],
        Data_CSV_sources: List[str]
) -> None:
    """
    Validate that all input lists have the same length.

    Raises:
        ValueError: If input lists have inconsistent lengths
    """
    lengths = {
        'chemNames': len(chemNames),
        'comptNames': len(comptNames),
        'ssType': len(ssType),
        'ssUnits': len(ssUnits),
        'Data_CSV_sources': len(Data_CSV_sources)
    }

    if len(set(lengths.values())) != 1:
        error_msg = '\n'.join(f'  {k}: {v}' for k, v in lengths.items())
        raise ValueError(
            f'Input data is not consistent. Number of elements:\n{error_msg}'
        )


def load_csv_data(csv_path: str) -> np.ndarray:
    """
    Load CSV data file.

    Args:
        csv_path: Path to CSV file

    Returns:
        Numpy array containing the data

    Raises:
        FileNotFoundError: If CSV file doesn't exist
    """
    csv_file = Path(csv_path)
    if not csv_file.exists():
        raise FileNotFoundError(f'CSV file not found: {csv_path}')

    return np.genfromtxt(csv_path, delimiter=',', skip_header=1)


def create_ss_json_structure(
        METADATA_Comment: str,
        METADATA_Source: str,
        chemNames: List[str],
        comptNames: List[str],
        ssType: List[str],
        ssUnits: List[str],
        Data_CSV_sources: List[str]
) -> Dict[str, Any]:
    """
    Create the sink-source JSON structure from input parameters.

    Args:
        METADATA_Comment: Comment for metadata
        METADATA_Source: Source for metadata
        chemNames: List of chemical names
        comptNames: List of compartment names
        ssType: List of sink-source types
        ssUnits: List of units
        Data_CSV_sources: List of CSV file paths

    Returns:
        Dictionary containing the complete JSON structure
    """
    # Initialize JSON structure
    ss_json_dic = {
        "METADATA": {
            "Comment": METADATA_Comment,
            "Source": METADATA_Source
        }
    }

    # Process each sink-source
    for idx, csv_source in enumerate(Data_CSV_sources):
        source_id = str(idx + 1)

        # Add general info with correct JSON key names
        ss_json_dic[source_id] = {
            "Chemical_name": chemNames[idx],
            "Compartment_name": comptNames[idx],
            "Type": ssType[idx],
            "Units": ssUnits[idx],
            "Data_Format": "JSON",
            "Data": {}
        }

        # Load and add CSV data
        ss_data = load_csv_data(csv_source)
        num_rows = ss_data.shape[0]

        for row_idx in range(num_rows):
            row_data = ss_data[row_idx, :].tolist()
            ss_json_dic[source_id]["Data"][str(row_idx + 1)] = row_data

    return ss_json_dic


def GenSSjson_fromCSV_driver(
        new_SS_JSONfile: str,
        METADATA_Comment: str,
        METADATA_Source: str,
        chemNames: List[str],
        comptNames: List[str],
        ssType: List[str],
        ssUnits: List[str],
        Data_CSV_sources: List[str]
) -> None:
    """
    Generate a sink-source JSON file from CSV data sources.

    Args:
        new_SS_JSONfile: Path for output JSON file
        METADATA_Comment: Comment for metadata
        METADATA_Source: Source for metadata
        chemNames: List of chemical names
        comptNames: List of compartment names
        ssType: List of sink-source types
        ssUnits: List of units
        Data_CSV_sources: List of CSV file paths containing the data
    """
    try:
        # Validate input consistency
        validate_input_consistency(
            chemNames, comptNames, ssType, ssUnits, Data_CSV_sources
        )

        # Create JSON structure
        ss_json_dic = create_ss_json_structure(
            METADATA_Comment,
            METADATA_Source,
            chemNames,
            comptNames,
            ssType,
            ssUnits,
            Data_CSV_sources
        )

        # Save to JSON file with compact array formatting
        output_path = Path(new_SS_JSONfile)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as outfile:
            outfile.write(compact_json_dump(ss_json_dic))

        print(f'✓ Successfully created: {new_SS_JSONfile}')

    except (ValueError, FileNotFoundError) as e:
        print(f'✗ ERROR: {e}')
        raise
    except Exception as e:
        print(f'✗ Unexpected error: {e}')
        raise


# Main execution
if __name__ == '__main__':
    # Configuration
    new_SS_JSONfile = 'examples/OpenWQ_source_example.json'

    METADATA_Comment = 'Test'
    METADATA_Source = 'test'

    chemNames = ["NO3", "NH4"]
    comptNames = ["SWE", "SWE"]
    ssType = ["source", "source"]
    ssUnits = ["kg", "kg"]
    Data_CSV_sources = [
        "examples/NO3_fertilizer.csv",
        "examples/NH4_fertilizer.csv"
    ]

    # Generate JSON file
    GenSSjson_fromCSV_driver(
        new_SS_JSONfile=new_SS_JSONfile,
        METADATA_Comment=METADATA_Comment,
        METADATA_Source=METADATA_Source,
        chemNames=chemNames,
        comptNames=comptNames,
        ssType=ssType,
        ssUnits=ssUnits,
        Data_CSV_sources=Data_CSV_sources
    )