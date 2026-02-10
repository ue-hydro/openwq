#!/bin/bash
# ============================================================================
# OpenWQ Supporting Scripts - Virtual Environment Setup
# ============================================================================
# This script creates a Python virtual environment with all dependencies
# needed to run the OpenWQ supporting scripts.
#
# Usage:
#   ./setup_venv.sh              # Create venv in default location (.venv)
#   ./setup_venv.sh /path/to/env # Create venv in custom location
#
# After setup, activate with:
#   source .venv/bin/activate    # Linux/macOS
#   .venv\Scripts\activate       # Windows
# ============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default venv location
VENV_DIR="${1:-.venv}"

echo "=============================================="
echo "OpenWQ Supporting Scripts - Environment Setup"
echo "=============================================="
echo ""

# Check Python version
PYTHON_CMD=""
for cmd in python3.11 python3.10 python3.9 python3 python; do
    if command -v $cmd &> /dev/null; then
        version=$($cmd -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
        major=$(echo $version | cut -d. -f1)
        minor=$(echo $version | cut -d. -f2)
        if [ "$major" -ge 3 ] && [ "$minor" -ge 8 ]; then
            PYTHON_CMD=$cmd
            break
        fi
    fi
done

if [ -z "$PYTHON_CMD" ]; then
    echo -e "${RED}ERROR: Python 3.8+ is required but not found.${NC}"
    echo "Please install Python 3.8 or later."
    exit 1
fi

echo -e "Using Python: ${GREEN}$PYTHON_CMD${NC} ($($PYTHON_CMD --version))"
echo -e "Virtual environment: ${GREEN}$VENV_DIR${NC}"
echo ""

# Create virtual environment
if [ -d "$VENV_DIR" ]; then
    echo -e "${YELLOW}Virtual environment already exists at $VENV_DIR${NC}"
    read -p "Do you want to recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        rm -rf "$VENV_DIR"
    else
        echo "Keeping existing environment. Updating packages..."
    fi
fi

if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment..."
    $PYTHON_CMD -m venv "$VENV_DIR"
fi

# Activate virtual environment
echo "Activating virtual environment..."
source "$VENV_DIR/bin/activate"

# Upgrade pip
echo ""
echo "Upgrading pip..."
pip install --upgrade pip

# Install requirements
echo ""
echo "Installing dependencies..."
pip install -r "$SCRIPT_DIR/requirements.txt"

# Verify installation
echo ""
echo "Verifying installation..."
python -c "
import sys
print(f'Python: {sys.version}')
print()
print('Checking dependencies:')

packages = [
    ('numpy', 'numpy'),
    ('pandas', 'pandas'),
    ('scipy', 'scipy'),
    ('h5py', 'h5py'),
    ('netCDF4', 'netCDF4'),
    ('geopandas', 'geopandas'),
    ('shapely', 'shapely'),
    ('fiona', 'fiona'),
    ('matplotlib', 'matplotlib'),
    ('requests', 'requests'),
]

all_ok = True
for name, module in packages:
    try:
        m = __import__(module)
        version = getattr(m, '__version__', 'unknown')
        print(f'  [OK] {name}: {version}')
    except ImportError as e:
        print(f'  [FAIL] {name}: {e}')
        all_ok = False

# Optional packages
print()
print('Optional packages:')
optional = [
    ('SALib', 'SALib'),
    ('jsonschema', 'jsonschema'),
]
for name, module in optional:
    try:
        m = __import__(module)
        version = getattr(m, '__version__', 'unknown')
        print(f'  [OK] {name}: {version}')
    except ImportError:
        print(f'  [SKIP] {name}: not installed')

if not all_ok:
    sys.exit(1)
"

if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}=============================================="
    echo "Setup complete!"
    echo "==============================================${NC}"
    echo ""
    echo "To activate the environment, run:"
    echo ""
    echo -e "  ${YELLOW}source $VENV_DIR/bin/activate${NC}"
    echo ""
    echo "To deactivate, run:"
    echo ""
    echo -e "  ${YELLOW}deactivate${NC}"
    echo ""
else
    echo ""
    echo -e "${RED}Setup completed with errors. Check the output above.${NC}"
    exit 1
fi
