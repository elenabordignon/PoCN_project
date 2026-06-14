#!/bin/bash
# run_grid_simulations_t4200.sh
# Script to run simulations starting at t=4200 in parallel for amplitudes:
# 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0

set -e

# Base directories
SIM_DIR="simulator"
DATA_DIR="data"

echo "=== STARTING PARALLEL GRID SIMULATIONS (t=4200) ==="

# Navigate to simulation directory
cd "$SIM_DIR"

# Force recompile of C++ simulator to pick up T_save = 4200 change
echo "Recompiling C++ simulator..."
g++ -std=c++11 -O3 -o tangled_nature_class tangled_nature_class.cpp

# Define the grid of amplitudes
AMPLITUDES=(0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
SEMENTE=46
FINESTRA=200

# Start all simulations in parallel
for amp in "${AMPLITUDES[@]}"; do
    echo "Starting background simulation for Amplitude: $amp (t=4200)"
    CARTELLA_OUT="../$DATA_DIR/raw_amp${amp}_t4200/"
    mkdir -p "$CARTELLA_OUT"
    
    # Run in the background using '&'
    ./tangled_nature_class $SEMENTE "$CARTELLA_OUT" $amp $FINESTRA > /dev/null 2>&1 &
done

echo "Waiting for all background simulations to complete..."
wait

echo "=== ALL PARALLEL GRID SIMULATIONS (t=4200) COMPLETED SUCCESSFULLY ==="
