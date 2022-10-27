#!/usr/bin/env bash

set -e

if ! command -v conda &> /dev/null
then
    echo ""
    echo "Conda could not be located on your system. Please follow installation"
    echo "instructions at https://docs.conda.io/en/latest/miniconda.html"
    echo "before running this script."
    echo ""
    exit
fi

# activate conda environment
eval "$(conda shell.bash hook)"

if conda env list | grep -q speedy_physics_compare; then
    # if speedy_physics_compare env already exists update it
    conda env update --file environment.yml --prune
else
    # create a new env
    conda env create -f environment.yml
fi

# Julia deps
conda activate speedy_physics_compare
# install julia dependecies
julia -e 'import Pkg; Pkg.add("PyCall")'
julia -e 'import Pkg; Pkg.develop(path="models/speedy_physics.jl")'
