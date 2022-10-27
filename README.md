# Speedy Physics Compare

This repository is used to evaluate ported [speedy_physics.f90](https://github.com/samhatfield/speedy.f90) to [SpeedyWeather.jl](https://github.com/milankl/SpeedyWeather.jl) parametrization functions.


## Requirements

- Linux/macOS/WSL
- Anaconda/Miniconda


## Set-up

Clone this repository with `--recurse-submodules` (i.e. `git clone --recurse-submodules https://github.com/dmey/speedy_physics_compare`) then, from your command line interface, type:

```sh
scripts/conda_init.sh
```


## Evaluation

You can view the current evaluation by opening [evaluate.ipynb](evaluate.ipynb). To run the evaluation, from your command line interface, type:

```sh
conda activate speedy_physics_compare
jupyter lab evaluate.ipynb
```
