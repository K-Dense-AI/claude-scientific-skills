#!/bin/bash
# QTL Analysis Dependencies Installer
#
# One-shot installer for all tools required by the qtl-analysis skill

set -e

echo "========================================"
echo "QTL Analysis Dependencies Installer"
echo "========================================"
echo ""

# Detect OS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS="linux"
    if command -v apt-get &> /dev/null; then
        PKG_MANAGER="apt"
    elif command -v yum &> /dev/null; then
        PKG_MANAGER="yum"
    elif command -v dnf &> /dev/null; then
        PKG_MANAGER="dnf"
    else
        PKG_MANAGER="unknown"
    fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macos"
    PKG_MANAGER="brew"
else
    echo "⚠️ Unknown OS: $OSTYPE"
    OS="unknown"
    PKG_MANAGER="unknown"
fi

echo "OS: $OS"
echo "Package Manager: $PKG_MANAGER"
echo ""

# Install Python packages
echo "Installing Python packages..."
pip install --upgrade pip

# Core dependencies
pip install pandas numpy matplotlib scipy scikit-learn statsmodels

# PyTorch (CPU version by default)
echo "Installing PyTorch (CPU)..."
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu

# Optional: GPU version
read -p "Install GPU support for PyTorch? (requires CUDA, y/n) [n]: " install_gpu
default="n"
install_gpu=${install_gpu:-$default}

if [[ $install_gpu == "y" ]]; then
    pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121
fi

# QTL-specific Python packages
echo "Installing tensorQTL and pyQTL..."
pip install tensorqtl qtl

echo "✅ Python packages installed"
echo ""

# Install CLI tools
if command -v conda &> /dev/null; then
    echo "Installing CLI tools via conda..."
    conda install -c bioconda -c conda-forge plink plink2 gemma r-qtl2 r-tidyverse -y
else
    echo "⚠️ Conda not found. Will attempt manual installation..."
    
    # PLINK
    if ! command -v plink &> /dev/null; then
        echo "Installing PLINK..."
        if [[ "$OS" == "linux" ]]; then
            wget https://s3.amazonaws.com/plink1-assets/dev/plink_linux_x86_64.zip
            unzip plink_linux_x86_64.zip -d ~/bin/
            rm plink_linux_x86_64.zip
        elif [[ "$OS" == "macos" ]]; then
            wget https://s3.amazonaws.com/plink1-assets/dev/plink_mac_x86_64.zip
            unzip plink_mac_x86_64.zip -d ~/bin/
            rm plink_mac_x86_64.zip
        fi
    fi
    
    # PLINK 2
    if ! command -v plink2 &> /dev/null; then
        echo "Installing PLINK 2..."
        if [[ "$OS" == "linux" ]]; then
            wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64.zip
            unzip plink2_linux_x86_64.zip -d ~/bin/
            rm plink2_linux_x86_64.zip
        elif [[ "$OS" == "macos" ]]; then
            wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_mac_x86_64.zip
            unzip plink2_mac_x86_64.zip -d ~/bin/
            rm plink2_mac_x86_64.zip
        fi
    fi
    
    # GEMMA
    if ! command -v gemma &> /dev/null; then
        echo "Installing GEMMA..."
        echo "Please install GEMMA manually from: https://github.com/genetics-statistics/GEMMA/releases"
    fi
    
    # R and qtl2
    if ! command -v R &> /dev/null; then
        echo "Installing R..."
        if [[ "$PKG_MANAGER" == "apt" ]]; then
            sudo apt-get update
            sudo apt-get install -y r-base
        elif [[ "$PKG_MANAGER" == "brew" ]]; then
            brew install r
        fi
    fi
    
    # R/qtl2
    echo "Installing R/qtl2..."
    Rscript -e 'install.packages("qtl2", repos="https://cloud.r-project.org")'
fi

echo ""
echo "========================================"
echo "Installation Complete"
echo "========================================"
echo ""
echo "Verify installation:"
echo "  python scripts/check_system.py"
echo ""
echo "Run examples:"
echo "  python examples/gwas-lmm/run_gwas.py"
echo "  python examples/eqtl-cis/run_eqtl.py"
echo "  python examples/classical-qtl/run_lodscan.py"
echo "  python examples/population-structure/run_structure.py"
