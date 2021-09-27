# A robust and efficient line search for self-consistent field iterations
[![](https://img.shields.io/badge/arxiv-2109.14018-red)](https://arxiv.org/abs/2109.14018)

Supporting information containing structures,
raw data and computational scripts for the paper:

Michael F. Herbst and Antoine Levitt  
*A robust and efficient line search for self-consistent field iterations*  
Preprint on [arxiv (2109.14018)](https://arxiv.org/abs/2109.14018)

The code in this repository has been used to run all calculations
and produce all plots of the above paper.
It relies on [DFTK](https://dftk.org) version 0.3.10
and [ASE](https://wiki.fysik.dtu.dk/ase/).

In particular:
  - The subfolders contain the employed structures and computational parameters
    of the numerical tests conducted for the paper. See in particular the
    `<case>/config.jl` file for the computational parameters of each test case
    (kinetic energy cutoff, k-Point sampling, smearing).
  - [analyse.jl](analyse.jl) contains the code to produce the plots of the paper.
  - [common.jl](common.jl) contains the routines to setup and run the calculations.
  - [plot_all_summaries.jl](plot_all_summaries.jl) plots the `summary.svg` plots
    of each subfolder which provide an overview of the convergence over all damping
    strategies considered.


## Running the code and reproducing the plots
Running the code requires an installation of
[Julia 1.6.0](https://julialang.org/downloads/#current_stable_release),
of [DFTK](https://docs.dftk.org/dev/guide/installation/)
and of [ASE](https://wiki.fysik.dtu.dk/ase/).
With this setup can generate the plots of the paper by executing:
```bash
julia --project=@. -e "import Pkg; Pkg.instantiate()"  # Install dependencies
julia run.jl      # Generate data
julia analyse.jl  # Generate plots
```

Be aware that generating the data takes a long time
(like two weeks on a couple of cluster nodes). The raw data is,
however, included in the repository, such that the plotting works
without running the calculations beforehand.

## List of test systems
  | Folder | System | Preconditioner | Anderson? | Comments |
  | ------ | ------ | -------------- | --------- | -------- |
  | [Al_nodiis](Al_nodiis)                 | Al8 supercell      | ×      | ×        | |
  | [Al_nodiis_Kerker](Al_nodiis_Kerker)   | Al8 supercell      | Kerker | ×        | |
  | [Al](Al)                               | Al40 supercell     | ×      | Anderson | |
  | [Al_Kerker](Al_Kerker)                 | Al40 supercell     | Kerker | Anderson | |
  | [AlVac](AlVac)                         | Al40 surface       | ×      | Anderson | |
  | ----                                   | ----               | ----   | ----     | ----           |
  | [GaAs](GaAs)                           | Ga40As40 supercell | ×      | Anderson | |
  | ----                                   | ----               | ----   | ----     | ----           |
  | [CoFeMnGa](CoFeMnGa)                   | CoFeMnGa           | Kerker | Anderson | |
  | [Fe2CrGa](Fe2CrGa)                     | Fe2CrGa            | Kerker | Anderson | |
  | [Fe2CrGa_nodiis](Fe2CrGa_nodiis)       | Fe2CrGa            | Kerker | ×        | |
  | [Fe2MnAl](Fe2MnAl)                     | Fe2MnAl            | Kerker | Anderson | |
  | [FeNiF6](FeNiF6)                       | FeNiF6             | Kerker | Anderson | |
  | [Mn2RuGa](Mn2RuGa)                     | Mn2RuGa            | Kerker | Anderson | |
  | [Mn3Si](Mn3Si)                         | Mn3Si              | Kerker | Anderson | |
  | ----                                   | ----               | ----   | ----     | ----           |
  | [Fe2CrGa_40guess](Fe2CrGa_40guess)     | Fe2CrGa            | Kerker | Anderson | started within the stagnating region |
  | [Fe2CrGa_maxcond02](Fe2CrGa_maxcond02) | Fe2CrGa            | Kerker | Anderson | using maximal conditioning 10² in Anderson |
  | [Fe2CrGa_maxcond04](Fe2CrGa_maxcond04) | Fe2CrGa            | Kerker | Anderson | using maximal conditioning 10⁴ in Anderson |
  | [Fe2CrGa_maxcond06](Fe2CrGa_maxcond06) | Fe2CrGa            | Kerker | Anderson | using maximal conditioning 10⁶ in Anderson |
  | [Fe2CrGa_maxcond08](Fe2CrGa_maxcond08) | Fe2CrGa            | Kerker | Anderson | using maximal conditioning 10⁸ in Anderson |
  | [Fe2CrGa_maxcond10](Fe2CrGa_maxcond10) | Fe2CrGa            | Kerker | Anderson | using maximal conditioning 10¹⁰ in Anderson |
