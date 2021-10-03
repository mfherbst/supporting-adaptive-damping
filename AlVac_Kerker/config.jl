lattice = load_lattice(joinpath(@__DIR__, "AlVac_2020.json"))
atoms   = load_atoms(  joinpath(@__DIR__, "AlVac_2020.json"))
atoms   = attach_pseudos(atoms, Al="hgh/pbe/al-q3.hgh")

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
    mixing           = KerkerMixing(),
    cases            = standard_cases(),
)
