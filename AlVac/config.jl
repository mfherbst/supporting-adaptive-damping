lattice = load_lattice(joinpath(@__DIR__, "AlVac_2020.json"))
atoms   = load_atoms(  joinpath(@__DIR__, "AlVac_2020.json"))
atoms   = attach_pseudos(atoms, Al="hgh/pbe/al-q3.hgh")

cases = standard_cases()
cases["adaptive_01"] = (damping=DFTK.AdaptiveDamping(0.1), )
cases["adaptive_05"] = (damping=DFTK.AdaptiveDamping(;α_trial_min=0.5, α_min=0.05), )

config[@__DIR__] = (
    directory        =  @__DIR__,
    lattice          = lattice,
    atoms            = atoms,
    magnetic_moments = [],
    #
    kwargs_model     = (smearing=Smearing.Gaussian(), temperature=0.001),
    Ecut             = 20,
    kgrid            = [6, 6, 1],
    supersampling    = 1.5,
    potential        = nothing,
    kwargs_scf       = (; ),
    mixing           = SimpleMixing(),
    cases            = cases,
)

