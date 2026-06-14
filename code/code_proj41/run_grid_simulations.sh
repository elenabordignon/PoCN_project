#!/bin/bash
# run_grid_simulations.sh
# Script to run simulations for the grid of perturbation amplitudes:
# 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0

set -e

# Base directories
BASE_DIR="/Users/elenabordignon/Desktop/CompNet/CN_exam/Tangled_nature_project"
SIM_DIR="$BASE_DIR/simulator"
DATA_DIR="$BASE_DIR/data"

echo "=== STARTING GRID SIMULATIONS ==="

# Navigate to simulation directory
cd "$SIM_DIR"

# Ensure the executable exists
if [ ! -f "./tangled_nature_class" ]; then
    echo "Compiling C++ simulator..."
    g++ -std=c++11 -O3 -o tangled_nature_class tangled_nature_class.cpp
fi

# Define the grid of amplitudes
AMPLITUDES=(0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
SEMENTE=46
FINESTRA=200

for amp in "${AMPLITUDES[@]}"; do
    echo "--------------------------------------------------"
    echo "Running simulation for Amplitude: $amp"
    CARTELLA_OUT="$DATA_DIR/raw_amp${amp}/"
    mkdir -p "$CARTELLA_OUT"
    
    # Run simulation
    ./tangled_nature_class $SEMENTE "$CARTELLA_OUT" $amp $FINESTRA
done

echo "=== ALL GRID SIMULATIONS COMPLETED SUCCESSFULLY ==="
