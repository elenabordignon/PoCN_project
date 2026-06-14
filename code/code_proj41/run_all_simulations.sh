#!/bin/bash
# run_all_simulations.sh
# Script to run simulations for missing perturbation amplitudes: 0.0, 0.25, 0.75

set -e

# Base directories
BASE_DIR="/Users/elenabordignon/Desktop/CompNet/CN_exam/Tangled_nature_project"
SIM_DIR="$BASE_DIR/simulator"
DATA_DIR="$BASE_DIR/data"

echo "=== STARTING SIMULATIONS ==="

# Make sure output directories exist
mkdir -p "$DATA_DIR/raw_amp0.0"
mkdir -p "$DATA_DIR/raw_amp0.25"
mkdir -p "$DATA_DIR/raw_amp0.75"

# Navigate to simulation directory
cd "$SIM_DIR"

# Ensure the executable exists
if [ ! -f "./tangled_nature_class" ]; then
    echo "Compiling C++ simulator..."
    g++ -std=c++11 -O3 -o tangled_nature_class tangled_nature_class.cpp
fi

# Run Amplitude 0.0 (no perturbation)
echo "Running simulation for Amplitude 0.0..."
./tangled_nature_class 46 "$DATA_DIR/raw_amp0.0/" 0.0 200

# Run Amplitude 0.25
echo "Running simulation for Amplitude 0.25..."
./tangled_nature_class 46 "$DATA_DIR/raw_amp0.25/" 0.25 200

# Run Amplitude 0.75
echo "Running simulation for Amplitude 0.75..."
./tangled_nature_class 46 "$DATA_DIR/raw_amp0.75/" 0.75 200

echo "=== ALL SIMULATIONS COMPLETED SUCCESSFULLY ==="
