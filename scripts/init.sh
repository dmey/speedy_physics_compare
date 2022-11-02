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
julia -e 'import Pkg; Pkg.add("PyCall")' # For PyJulia
julia -e 'import Pkg; Pkg.add("DataStructures")' # For fortran2julia.jl
julia -e 'import Pkg; Pkg.develop(path="models/speedy_physics.jl")'


# Generate python bindings for fortran funcs

# Shortwave radiation
pushd models/speedy_physics.f90/src
    rm -f *.so
    MODULE_NAME=radiation.f90
    MODULE_NAME_TMP=radiation.temp.f90
    cp -f $MODULE_NAME $MODULE_NAME_TMP
    # Some funcs are not exported
    sed -i '/public sol_oz, cloud, radsw/c\public sol_oz, cloud, radsw, solar, solar_diurnal' $MODULE_NAME_TMP
    f2py -c -m shortwave_radiation $MODULE_NAME_TMP
    rm -f $MODULE_NAME_TMP
popd

# Longwave radiation
pushd models/speedy_physics.f77
    rm -f *.so
    
    MODULE_NAME=phy_radiat.f
    MODULE_NAME_TMP=phy_radiat.temp.f

    # Downwelling: patch intent and create bindings
    cp -f $MODULE_NAME $MODULE_NAME_TMP
    sed -i '/NL1=NLEV-1/i \
      INTEGER IMODE\
Cf2py intent(in) IMODE, TA, TS\
Cf2py intent(out) FSFCD, FSFCU, DFABS\
Cf2py intent(out) FSFC, FTOP\
      call radset' $MODULE_NAME_TMP
    f2py -c -m longwave_radiation_down $MODULE_NAME_TMP

    # Upwelling: patch intent and create bindings
    cp -f $MODULE_NAME $MODULE_NAME_TMP
    sed -i '/NL1=NLEV-1/i \
      INTEGER IMODE\
Cf2py intent(in) IMODE, TA, TS, FSFCD, FSFCU, DFABS\
Cf2py intent(out) FSFC, FTOP\
      call radset' $MODULE_NAME_TMP
    f2py -c -m longwave_radiation_up $MODULE_NAME_TMP

    rm -f $MODULE_NAME_TMP
popd
