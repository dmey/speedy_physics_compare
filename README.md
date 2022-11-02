# Speedy Physics Compare

This repository is used to evaluate ported [speedy_physics.f90](https://github.com/samhatfield/speedy.f90) to [SpeedyWeather.jl](https://github.com/milankl/SpeedyWeather.jl) parametrization functions.


## Requirements

- Linux/macOS/WSL
- [Anaconda/Miniconda](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links)


## Set-up

Clone this repository with `--recurse-submodules` (i.e. `git clone --recurse-submodules https://github.com/dmey/speedy_physics_compare`) then, to install the required dependencies and generate python bindings for fortran modules type the following command from your command line interface:

```sh
scripts/init.sh
```


## Evaluation

You can view the current evaluation by opening the evaluation notebooks (`evaluate_`). Alternatively, to run the evaluation locally type the following commands from your command line interface, :

```sh
conda activate speedy_physics_compare
jupyter lab evaluate.ipynb
```
