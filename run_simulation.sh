#!/bin/bash

# Navigate to the project root (assuming the script is in the root)
# If you run this from elsewhere, adjust the path or comment this out.
# cd "$(dirname "$0")"

# Ensure the executable exists
if [ ! -f bin/simulation ]; then
    echo "Executable not found. Please run 'make' first from the project root."
    exit 1
fi

# Run the simulation with example arguments
# Replace with your actual input files and desired parameters
# Example parameters based on your main function logic:
bin/simulation \
    -e 1 \
    -d 1 \
    -g 1 \
    -t 1 \
    -i data/your_input_sigma_file.dat \
    -s 0.01 \
    # Add any other parameters your parse_options function expects (e.g., -n for NGRID if optinp == 1)

echo "Simulation finished."
